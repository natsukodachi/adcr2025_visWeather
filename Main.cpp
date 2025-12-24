# include <Siv3D.hpp>
# include <netcdf>

using namespace netCDF;

// ERA5 の PMSL（海面更正気圧）を格納する構造体
struct weatherData
{
	Grid<double> pmsl;  // 海面更正気圧 [hPa]
	Array<double> lats; // 緯度配列 [度]
	Array<double> lons; // 経度配列 [度]
};

// NetCDF から PMSL データを読み込む
void LoadPMSL(weatherData& field, const FilePath& ncPath)
{
	// NetCDF ファイルを読み込みモードで開く
	NcFile nc(Unicode::Narrow(ncPath), NcFile::read);

	// 変数ハンドルを取得
	NcVar latVar = nc.getVar("latitude");
	NcVar lonVar = nc.getVar("longitude");
	NcVar mslVar = nc.getVar("msl"); // ERA5 の海面更正気圧（単位: Pa）

	// 必須変数が存在するかチェック
	if (latVar.isNull() || lonVar.isNull() || mslVar.isNull())
	{
		throw Error{ U"latitude / longitude / msl が見つかりません" };
	}

	// 次元サイズ取得（緯度・経度）
	const size_t nLat = latVar.getDim(0).getSize();
	const size_t nLon = lonVar.getDim(0).getSize();


	// 緯度・経度配列を読み込み
	field.lats.resize(nLat);
	field.lons.resize(nLon);
	latVar.getVar(field.lats.data());
	lonVar.getVar(field.lons.data());

	// msl の読み込み範囲指定（time/level/lat/lon の順を想定）
	Array<size_t> start = { 0, 0, 0 };
	Array<size_t> count = { 1, nLat, nLon };

	// バッファに読み込み（単位: Pa）
	Array<float> buf(nLat * nLon);
	mslVar.getVar(start, count, buf.data());

	// Grid は pmsl[y][x] なので width=nLon, height=nLat に整形
	field.pmsl.resize(nLon, nLat); // width=nLon, height=nLat（Grid は pmsl[y][x]）

	// Pa → hPa に換算して格納
	for (size_t j = 0; j < nLat; ++j)
	{
		for (size_t i = 0; i < nLon; ++i)
		{
			const double pa = buf[j * nLon + i];
			field.pmsl[j][i] = pa * 0.01; // 1 hPa = 100 Pa
		}
	}
}

// カラーマップ画像を作成（pmsl 用）
Image CreateColormapImage(const Grid<double>& data, double vmin, double vmax, ColormapType cmapType)
{
	// 値の範囲を 0.0-1.0 に正規化
	const double invRange = 1.0 / (vmax - vmin);
	const int32 w = static_cast<int32>(data.width());
	const int32 h = static_cast<int32>(data.height());

	Image img(w, h);
	for (int32 y = 0; y < h; ++y)
	{
		for (int32 x = 0; x < w; ++x)
		{
			double t = (data[y][x] - vmin) * invRange;
			t = Clamp(t, 0.0, 1.0);
			img[y][x] = Colormap01F(t, cmapType);
		}
	}
	return img;
}

// 表示範囲（経度・緯度の min/max）
struct GeoBounds
{
	double lonMin, lonMax;
	double latMin, latMax;
};

// 配列から min/max を求める
static GeoBounds GetBounds(const weatherData& field)
{
	GeoBounds bounds{};
	bounds.lonMin = Math::Inf;
	bounds.lonMax = -Math::Inf;
	for (const double lon : field.lons)
	{
		bounds.lonMin = Min(bounds.lonMin, lon);
		bounds.lonMax = Max(bounds.lonMax, lon);
	}

	bounds.latMin = Math::Inf;
	bounds.latMax = -Math::Inf;
	for (const double lat : field.lats)
	{
		bounds.latMin = Min(bounds.latMin, lat);
		bounds.latMax = Max(bounds.latMax, lat);
	}
	return bounds;
}

// Grid の最小値・最大値を取得（min は 100 以上の値のみを対象）
std::pair<double, double> GetMinMax(const Grid<double>& data) {
	double mn = std::numeric_limits<double>::infinity();
	double mx = -std::numeric_limits<double>::infinity();
	bool hasMin = false; // 100 以上の有効値が見つかったか
	const int32 w = static_cast<int32>(data.width());
	const int32 h = static_cast<int32>(data.height());
	for (int32 y = 0; y < h; ++y) {
		for (int32 x = 0; x < w; ++x) {
			const double v = data[y][x];
			// 最小値は 100 以上の値のみを考慮
			if (v >= 100.0) {
				if (!hasMin || v < mn) mn = v;
				hasMin = true;
			}
			// 最大値は従来通り全値を対象
			if (v > mx) mx = v;
		}
	}
	// 100 以上が見つからなかった場合のフォールバック
	if (!hasMin) {
		mn = 100.0;
	}
	return { mn, mx };
}

// 海岸線
class CoastlineOverlay
{
public:
	// 可視範囲を受け取り、範囲に入る国だけ抽出
	CoastlineOverlay(const RectF& geoViewRectLonY)
	{

		countries_ = GeoJSONFeatureCollection{ JSON::Load(U"example/geojson/countries.geojson") }
			.getFeatures()
			.map([](const GeoJSONFeature& f) { return f.getGeometry().getPolygons(); });

		// 可視範囲内の国だけ抽出
		visibleIndices_.reserve(countries_.size());
		for (auto&& [i, country] : Indexed(countries_))
		{
			if (country.computeBoundingRect().intersects(geoViewRectLonY))
			{
				visibleIndices_ << i;
			}
		}
	}

	// テクスチャの描画領域（destRect）に合わせて海岸線を重ね描き
	void draw(const RectF& destRect, double lonMin, double lonMax, double yMin, double yMax) const
	{
		// グリッドの 1 ピクセルサイズ
		const Vec2 pixelSize = destRect.size / Vec2{ 1.0 * m_w, 1.0 * m_h };
		const Vec2 halfPixel = pixelSize * 0.5;

		// 経度・緯度（y は-latitude）→ピクセル座標へのスケール
		const double sx = (destRect.w - pixelSize.x) / (lonMax - lonMin); // (w-1) 相当
		const double sy = (destRect.h - pixelSize.y) / (yMax - yMin);     // (h-1) 相当

		// 左上原点合わせ（半ピクセル補正込み）
		const Vec2 translate = (destRect.pos + halfPixel) - Vec2{ lonMin * sx, yMin * sy };

		// 座標変換の適用
		const Transformer2D t{ Mat3x2::Scale(Vec2{ sx, sy }).translated(translate) };

		// 画面拡大率に応じた線幅（スケール非依存）
		const double lineThickness = (1.0 / Graphics2D::GetMaxScaling());

		// 可視国のみ描画
		for (const size_t i : visibleIndices_)
		{
			countries_[i].drawFrame(lineThickness, ColorF{ 0.1, 0.1, 0.1, 0.85 });
		}
	}

	// グリッドサイズ（ピクセル数）を設定
	void setGridSize(int32 w, int32 h)
	{
		m_w = w;
		m_h = h;
	}

private:
	Array<MultiPolygon> countries_;   // 国境線ポリゴン
	Array<size_t> visibleIndices_;    // 可視範囲に入る国のインデックス
	int32 m_w = 1;                    // グリッド幅
	int32 m_h = 1;                    // グリッド高さ
};

void Main()
{
	// ウィンドウ初期化
	Window::Resize(600, 600);
	Scene::SetBackground(ColorF{ 0.2, 0.3, 0.4 });

	// NetCDF 入力ファイル（例）
	const String ncPath = U"pmsl.nc";
	weatherData field;
	LoadPMSL(field, ncPath); // PMSL を読み込み

	const int32 w = static_cast<int32>(field.pmsl.width());
	const int32 h = static_cast<int32>(field.pmsl.height());

	const auto [pmin, pmax] = GetMinMax(field.pmsl);
	Print << U"pmsl min:" << pmin << U" pmsl max:" << pmax;

	// カラーマップ用の表示範囲をデータの最小値・最大値に合わせる（単位: hPa）
	const double vmin = pmin;
	const double vmax = pmax;
	Texture cmapTex{ CreateColormapImage(field.pmsl, vmin, vmax, ColormapType::Turbo) };

	// 経度・緯度の範囲取得
	const auto bounds = GetBounds(field);

	// GeoJSON の座標系が y = -latitude 前提なので変換
	const double yMin = -bounds.latMax;
	const double yMax = -bounds.latMin;

	// 可視範囲矩形（経度 x、-緯度 y）
	const RectF geoViewRect{ bounds.lonMin, yMin, (bounds.lonMax - bounds.lonMin), (yMax - yMin) };

	// 海岸線オーバーレイの準備
	CoastlineOverlay coastline{ geoViewRect };
	coastline.setGridSize(w, h);

	// メインループ
	while (System::Update())
	{
		// テクスチャをウィンドウに収まるようにスケール
		const double drawScale = Min(Scene::Width() / static_cast<double>(w), Scene::Height() / static_cast<double>(h));
		cmapTex.scaled(drawScale).drawAt(Scene::CenterF()); // 中央に描画

		// 描画先矩形（テクスチャ領域）
		const Vec2 drawSize = Vec2{ w, h } * drawScale;
		const RectF destRect{ Scene::CenterF() - (drawSize * 0.5), drawSize };

		// 海岸線を重ね描き
		coastline.draw(destRect, bounds.lonMin, bounds.lonMax, yMin, yMax);
	}
}
