from pathlib import Path

from postprocess.utils_plotter import create_folder
from postprocess.visualize import Visualize

main = Path.cwd()
push_dir = main / "RCMRF/Medium"
create_folder(push_dir / "figs")

f = push_dir / "SPO_x.json"
spo_model = f

# # Labels if multiple graphs in one figure
# labels = ["SPO"]
# # For figure naming
# name = "x_medium"

# viz = Visualize(export=False, filetype="pdf", export_dir=push_dir / "figs", flag=True)
# viz.plot_spo(spo_model, name=name, labels=labels)

filename = [main / "RCMRF/Medium/SPO_x.json",
            main / "RCMRF/Medium/SPO_y.json",
            main / "RCMRF/High/SPO_x.json",
            main / "RCMRF/High/SPO_y.json"]
name = "all"
labels = ["Medium, x", "Medium, y", "High, x", "High, y"]

viz = Visualize(export=True, filetype="pdf", export_dir=push_dir / "figs", flag=True)
viz.plot_spo(filename, name=name, labels=labels)
