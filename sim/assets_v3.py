"""Pitch-deck assets: modulation-vs-depth chart (PNG + SVG) and a decoder demo GIF.

Both are generated from simulation output; both carry a visible label saying so.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

import config_v3 as C  # noqa: E402
import decode_sim as DS  # noqa: E402
import decoder  # noqa: E402
import make_tables_v3 as M3  # noqa: E402
import run_v2 as R2  # noqa: E402
import timesim  # noqa: E402

BG, INK = "#F6F4EF", "#0B0F14"
LINES = {30: "#FF4D2E", 100: "#B23A6B", 300: "#5B6573"}
MUTED = "#5B6573"


def decoder_floor_60s():
    fl = M3.load_csv(os.path.join(M3.V2, "decoder_floor.csv"))
    v = [f["m_min"] for f in fl if f["reference"] and f["config"].startswith("scan") and f["read_s"] > 50]
    return max(v)


def chart(res, out):
    floor = decoder_floor_60s()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=200)
    fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
    depths = [0] + C.DEPTHS_V2_UM
    for a, col in LINES.items():
        m, e = [], []
        for d in depths:
            r = M3.one(res, placement="A_top" if d == 0 else "C_embedded", depth_um=d, spot_um=a, cam="scan_NA0.25",
                       wavelength_nm=635, musp_780=1.0, **M3.BASE)
            m.append(r["m"]); e.append(r.get("m_err", 0.0) if np.isfinite(r.get("m_err", np.nan)) else 0.0)
        ax.errorbar(depths, m, yerr=e, color=col, lw=2.2, marker="o", ms=4, capsize=2.5, elinewidth=1)
        ax.annotate(f"{a} µm SPOT", xy=(depths[-1], m[-1]), xytext=(8, 0), textcoords="offset points",
                    color=col, va="center", fontsize=11, fontweight="bold")
    ax.axhline(floor, color=INK, lw=1, ls=(0, (4, 3)))
    ax.text(5, floor * 1.25, f"decoder floor, 60 s read (m ≈ {floor:.1e}); chips deeper than 500 µm need ~5× more",
            color=INK, fontsize=9)
    ax.set_yscale("log")
    ax.set_xlim(-20, 1000 + 20); ax.set_ylim(1e-5, 0.6)
    ax.set_xlabel("chip depth below the organoid surface (µm)", color=INK)
    ax.set_ylabel("modulation depth (blink contrast)", color=INK)
    ax.set_title("How much of the light blinks, by chip depth", color=INK, loc="left", fontsize=13, fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK)
    ax.tick_params(colors=INK)
    ax.grid(axis="y", color=INK, alpha=0.08, which="major")
    fig.text(0.01, 0.01, "Simulation, not measurement. 635 nm (measured reflectance change), average tissue, "
             "NA 0.25 scanning reader, 10 mW/cm².", color=MUTED, fontsize=8)
    fig.subplots_adjust(left=0.11, right=0.83, top=0.9, bottom=0.17)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"chart_modulation_vs_depth.{ext}"), facecolor=BG)
    plt.close(fig)


def _font(size, bold=False):
    path = os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf")
    return ImageFont.truetype(path, size)


def gif(res, out, seed=11):
    """Synthetic frames for 3 chips lying on the dish, read by a 20 fps plate camera (blink 5 Hz)."""
    row = M3.one(res, placement="B_under", depth_um=0, spot_um=100, cam="plate20MP_NA0.02", wavelength_nm=635,
                 musp_780=1.0, **M3.BASE)
    b = R2.budget_for(row)
    cam = b["cam"]
    rng = np.random.default_rng(seed)
    p = cam["obj_pixel_um"]
    pts = rng.normal(0, 0.45 * 100 / 2, (4000, 2))           # ~100 um SPOT, slightly blurred
    truth = [(17, 42), (99, 200), (250, 7)]
    sites = [(5.0, 5.0), (5.0, 15.0), (14.0, 10.0)]
    chips = [dict(pos_um=(x * p + p / 2, y * p + p / 2), pts_um=pts, e_mean=0.5 * (b["e_on"] + b["e_off"]),
                  e_mod=b["e_on"] - b["e_off"], id=i, value=v) for (i, v), (y, x) in zip(truth, sites)]
    n_c = 2
    sc = timesim.make_scene(chips, (20, 20), p, cam["fps"], cam["blink_hz"], n_c, b["e_bg_um2"], cam["read_noise"],
                            rng, sigma_src=1e-3)
    T = (C.FRAME_BITS + 2) * n_c / cam["blink_hz"]
    frames = timesim.render(sc, T, rng)
    dec = decoder.decode(frames, cam["fps"], cam["blink_hz"], n_c, use_reference=True, radius_px=3)
    lo, hi = np.percentile(frames, 1), np.percentile(frames, 99.7)
    S = 18                                   # screen pixels per camera pixel
    W, Hh = 20 * S + 300, 20 * S + 70
    f_big, f_small, f_bold = _font(15), _font(12), _font(15, True)
    imgs = []
    n = frames.shape[0]
    for k in range(n):
        g = np.clip((frames[k] - lo) / (hi - lo), 0, 1)
        g = (40 + 215 * g).astype(np.uint8)
        im = Image.new("RGB", (W, Hh), BG)
        cam_img = Image.fromarray(g, "L").resize((20 * S, 20 * S), Image.NEAREST).convert("RGB")
        im.paste(cam_img, (10, 50))
        d = ImageDraw.Draw(im)
        d.text((10, 10), "SYNTHETIC DATA", fill="#FF4D2E", font=f_bold)
        d.text((170, 12), f"simulated camera, 20 fps   t = {k / cam['fps']:5.1f} s", fill=INK, font=f_small)
        done = k >= n - 1
        prog = (k + 1) / n
        x0 = 20 * S + 25
        d.text((x0, 50), "decoder output", fill=INK, font=f_bold)
        for i, ((cid, val), (yy, xx)) in enumerate(zip(truth, sites)):
            cx, cy = 10 + int((xx + 0.5) * S), 50 + int((yy + 0.5) * S)
            d.ellipse([cx - 2.2 * S, cy - 2.2 * S, cx + 2.2 * S, cy + 2.2 * S], outline="#FF4D2E", width=2)
            d.text((cx + 2.4 * S, cy - 8), f"#{i + 1}", fill="#FF4D2E", font=f_small)
            y0 = 85 + 60 * i
            if done:
                m = [q for q in dec if abs(q["x"] - xx) < 2.5 and abs(q["y"] - yy) < 2.5]
                txt = f"chip {i + 1}:  ID {m[0]['id']}   value {m[0]['value']}" if m and m[0]["preamble_ok"] else f"chip {i + 1}: not decoded"
                d.text((x0, y0), txt, fill=INK, font=f_big)
            else:
                d.text((x0, y0), f"chip {i + 1}:  listening…", fill=MUTED, font=f_big)
                d.rectangle([x0, y0 + 24, x0 + 240, y0 + 32], outline=MUTED)
                d.rectangle([x0, y0 + 24, x0 + int(240 * prog), y0 + 32], fill="#B23A6B")
        d.text((10, Hh - 18), "Frames generated from simulated photon budgets + noise model; decoder never sees the truth.",
               fill=MUTED, font=f_small)
        imgs.append(im)
    imgs += [imgs[-1]] * 40                     # hold the decoded result ~2 s
    imgs[0].save(os.path.join(out, "decoder_demo.gif"), save_all=True, append_images=imgs[1:], duration=40, loop=0,
                 optimize=True)
    strip = Image.new("RGB", (W * 3, Hh), BG)
    for i, k in enumerate((0, n // 2, len(imgs) - 1)):
        strip.paste(imgs[k], (i * W, 0))
    strip.save(os.path.join(out, "decoder_demo_strip.png"))
    ok = sum(1 for (cid, val) in truth for q in dec if q["id"] == cid and q["value"] == val and q["preamble_ok"])
    return dict(decoded=ok, chips=len(truth), seconds=T)


def main(res, out):
    chart(res, out)
    g = gif(res, out)
    print(f"assets written; GIF decoded {g['decoded']}/{g['chips']} chips from {g['seconds']:.1f} s of synthetic frames")
    return g


if __name__ == "__main__":
    import pickle
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "v3")
    main(pickle.load(open(os.path.join(out, "mc_v3.pkl"), "rb"))["res"], out)
