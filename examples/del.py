from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from analysis.multiStripeAnalysis import get_records, get_ground_motion_batches, get_ground_motion
from utils.utils import read_text_file


def _append_record(x, y):
	"""
	Record selection issue related to the database of records
	Whereby the two pairs of the records have different sizes
	To remedy the issue, the smallest of the records is appended with zeros
	:return: list, list
	"""
	nx = len(x)
	ny = len(y)

	if nx < ny:
		x = np.append(x, np.zeros(ny - nx))
	if ny < nx:
		y = np.append(y, np.zeros(nx - ny))
	return x, y


seism = "Medium"
gmdir = Path.cwd().parents[2] / f"ProjectsCurrent/NSPerformanceAssessment/sample/groundMotionMSA/{seism}"
gmfileNames = ["GMR_H1_names.txt", "GMR_H2_names.txt", "GMR_dts.txt"]
outputsDir = Path.cwd() / "outputs/msa"

records = get_records(gmdir, outputsDir / "MSA", gmfileNames)

item = list(records.items())

for i in item:
	name, data = i
	print("Name", name)
	names_x = data["X"]
	names_y = data["Y"]
	dts = data["dt"]

	for rec in range(len(names_x)):
	# for rec in [24]:
		print(rec)
		eq_name_x = gmdir / name / names_x[rec]
		eq_name_y = gmdir / name / names_y[rec]
		dt = dts[rec]

		accg_x = read_text_file(eq_name_x)
		accg_y = read_text_file(eq_name_y)
		accg_x, accg_y = _append_record(accg_x, accg_y)

		dur = round(10 + dt * len(accg_x), 5)

		eq_x = accg_x
		eq_y = accg_y

		# Create time history
		# time_history = np.arange(dt, round(dur - 10, 5), dt)
		time_history = np.linspace(dt, round(dur - 10 + dt, 5), int(round((dur - 10) / dt, 0)))
		print(len(accg_x), len(accg_y), len(time_history), dt, (dur - 10) / dt)

		function_x = interp1d(time_history, eq_x, bounds_error=False, fill_value=0)
		function_y = interp1d(time_history, eq_y, bounds_error=False, fill_value=0)

		# try:
		# 	time_history = np.arange(dt, round(dur - 10 + dt, 5), dt)
		# 	print(len(accg_x), len(accg_y), len(time_history), dt, round(dur - 10 + dt, 5))
		# 	# Create interpolation function
		# 	function_x = interp1d(time_history, eq_x, bounds_error=False, fill_value=0)
		# 	function_y = interp1d(time_history, eq_y, bounds_error=False, fill_value=0)
		# except:
		# 	time_history = np.arange(dt, round(dur - 10, 5), dt)
		# 	print("t", len(accg_x), len(accg_y), len(time_history), dt, round(dur - 10, 5))
		# 	# Create interpolation function
		# 	function_x = interp1d(time_history, eq_x, bounds_error=False, fill_value=0)
		# 	function_y = interp1d(time_history, eq_y, bounds_error=False, fill_value=0)


exit()
time = np.arange(dt_gm, dur - dt_gm, dt_gm/4)
acc = function_x(time)

plt.plot(time, acc, color="b")
plt.plot(time_history, accg_x, color="r", ls="--")
plt.show()
