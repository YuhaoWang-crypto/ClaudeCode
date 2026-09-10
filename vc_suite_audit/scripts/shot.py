import sys
from playwright.sync_api import sync_playwright
CHROME="/opt/pw-browsers/chromium-1194/chrome-linux/chrome"
path=sys.argv[1]; out=sys.argv[2]; theme=sys.argv[3] if len(sys.argv)>3 else "light"
errs=[]
with sync_playwright() as p:
    b=p.chromium.launch(executable_path=CHROME, args=["--no-sandbox"])
    pg=b.new_page(viewport={"width":1180,"height":1500}, color_scheme=theme)
    pg.on("console", lambda m: errs.append(f"{m.type}: {m.text}") if m.type=="error" else None)
    pg.on("pageerror", lambda e: errs.append("PAGEERROR: "+str(e)))
    pg.goto("file://"+path); pg.wait_for_timeout(1800)
    print("tbody rows:", pg.eval_on_selector_all("table tbody tr","n=>n.length"))
    print("verdict cards:", pg.eval_on_selector_all(".vc","n=>n.length"))
    print("polylines:", pg.eval_on_selector_all("svg polyline","n=>n.length"))
    print("h-overflow:", pg.evaluate("document.documentElement.scrollWidth>document.documentElement.clientWidth"))
    pg.screenshot(path=out, full_page=True)
    b.close()
print("JS errors:", errs[:8] or "none")
