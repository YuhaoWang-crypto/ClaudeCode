"""CJK font provisioning for the Chinese-language report and figures.

The container ships WenQuanYi Zen Hei as a .ttc collection. matplotlib and
reportlab are both happier with a plain .ttf, so the first subfont is extracted
once into a cache directory. The extracted file is deliberately NOT committed
(it is ~11 MB and regenerable); callers get a path and the registered family
name.

WenQuanYi Zen Hei ships a Regular weight only, so the report builds its visual
hierarchy from size and colour rather than from bold.
"""

from __future__ import annotations

import pathlib
import shutil

TTC_CANDIDATES = [
    "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
    "/usr/share/fonts/truetype/arphic/uming.ttc",
]
TTF_CANDIDATES = [
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttf",
]
CACHE = pathlib.Path("assets/fonts/wqy-zenhei.ttf")


def ensure_cjk_ttf(cache: pathlib.Path | str = CACHE) -> pathlib.Path:
    """Return a path to a usable CJK .ttf, extracting from a .ttc if needed."""
    cache = pathlib.Path(cache)
    if cache.exists() and cache.stat().st_size > 0:
        return cache
    cache.parent.mkdir(parents=True, exist_ok=True)

    for ttf in TTF_CANDIDATES:
        p = pathlib.Path(ttf)
        if p.exists() and p.suffix == ".ttf":
            shutil.copy(p, cache)
            return cache

    for ttc in TTC_CANDIDATES:
        p = pathlib.Path(ttc)
        if not p.exists():
            continue
        from fontTools.ttLib import TTCollection
        coll = TTCollection(str(p))
        font = coll.fonts[0]
        font.flavor = None
        font.save(str(cache))
        return cache

    raise RuntimeError(
        "No CJK font found. Install one (e.g. fonts-wqy-zenhei) or place a "
        f".ttf at {cache}.")


def register_matplotlib(cache: pathlib.Path | str = CACHE) -> str:
    """Register the CJK font with matplotlib; returns the family name."""
    import matplotlib
    from matplotlib import font_manager
    path = ensure_cjk_ttf(cache)
    font_manager.fontManager.addfont(str(path))
    name = font_manager.FontProperties(fname=str(path)).get_name()
    matplotlib.rcParams["font.sans-serif"] = [name]
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["axes.unicode_minus"] = False
    return name


def register_reportlab(cache: pathlib.Path | str = CACHE,
                       name: str = "CJK") -> str:
    """Register the CJK font with reportlab; returns the font name.

    Bold/italic slots map back to the same face because the font has only one
    weight -- this keeps <b> markup from raising, it does not embolden.
    """
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    path = ensure_cjk_ttf(cache)
    pdfmetrics.registerFont(TTFont(name, str(path)))
    pdfmetrics.registerFontFamily(name, normal=name, bold=name,
                                  italic=name, boldItalic=name)
    return name
