# Screenshots for USER-GUIDE.md

The images are 600 dpi PNG captures of the plugin window, made with the plugin's own capture
(double-click the status bar, `widgets.save_widget_image`, default `dpi=600`).

To capture or refresh them:

1. Install the plugin and restart PyMOL. Open Plugin → pyPackmol.
2. Keep the default window size (820 × 680), so all images have the same size.
3. Use the example from section 10 of the guide (water and urea in a box), so the tabs show content:
   load `mixture.inp`, run it once in the Pack tab, click Get on `water` in the Sandbox, and so on.
4. For each image below, switch to the tab, double-click the status bar and save under the file name
   shown, in this folder.

| File | Tab | Suggested content |
|---|---|---|
| `setup-basic.png` | Set Up Input → Basic | `mixture.inp` loaded, a constraint rule selected, a volume check report. |
| `setup-additional.png` | Set Up Input → Additional | Any option selected, e.g. `maxit`. |
| `pack.png` | Pack | Log of a finished run, progress at 100%. |
| `sandbox.png` | Sandbox | A few rows from Get (Baricenter, Box, Sphere) on `water`. |
| `volume.png` | Volume Calculator | History with a molecule volume, a CGO volume and a `try` result. |
| `restart.png` | Restart Files | Log of a finished Partial load, ending with `EOF`. |
| `settings.png` | Plugin Settings | Default values. |

The guide shows each image 700 px wide (`<img … width="700">`), so the full resolution is only
visible when the image is opened on its own.
