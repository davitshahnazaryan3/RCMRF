from pathlib import Path

from postprocess.utils_plotter import create_folder
from postprocess.visualize import Visualize

main = Path.cwd()
push_dir = main / "outputs/pushover"
create_folder(push_dir / "figs")

f = push_dir / "SPO.json"
spo_model = f

# Labels if multiple graphs in one figure
labels = ["SPO"]
# For figure naming
name = "spo"

viz = Visualize(export=False, filetype="pdf", export_dir=push_dir / "figs", flag=True)
viz.plot_spo(spo_model, name=name, labels=labels)
