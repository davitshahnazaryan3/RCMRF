<h1 align="center">Generic RC MRF model creator</h1>

Procedures to create a reinforced concrete building (3D), or a frame (2D) nonlinear building model
consisting of moment resisting frames (MRF) as the primary lateral load resisting system.

May use multiprocessing to carry out multiple stripe analysis (non-linear time history analyses).

May be used in any os system.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5753463.svg)](https://doi.org/10.5281/zenodo.5753463)

**Required libraries**: requirements.txt

      python -m pip install -r requirements.txt

### Table of Contents
<details open>
<summary>Show/Hide</summary>
<br>

1. [Supported analysis types](#analysis)
2. [Required input arguments](#input)
3. [Samples](#samples)
4. [Future](#future)

</details>


### Supported analysis types
<details>
<a name="analysis"></a>
<summary>Show/Hide</summary>
<br>

1. Static elastic analysis (ST)
2. Modal analysis (MA)
3. Equivalent lateral force (ELF)
4. Static pushover analysis (PO)
5. Non-linear time history analysis (NLTHA)
   1. Incremental dynamic analysis (IDA)
	2. Multiple stripe analysis (MSA)
	
</details>


### Required input arguments per analysis
<details>
<a name="input"></a>
<summary>Show/Hide</summary>
<br>

* **sections_file** - required for All
  	
		csv or pickle file or DataFrame with hysteretic model parameters

* **loads_file** - required for All

		csv file masses and gravity loads

* **materials_file** - required for All

		csv file containing material properties

* **outputsDir** - required for All

		Directory to export outputs to

* **gmdir** - required for IDA and MSA
		
		Directory to read records from

* **gmfileNames** - required for IDA and MSA

		File names in order of ["GM_names_x", "GM_names_y", "GM_time_step"]

* **IM_type** - required for IDA, default to 2

		Intensity measure type

* **max_runs** - required for IDA, default to 15

		Maxium number of runs per record

* **analysis_time_step** - required for IDA and MSA, default to 0.01

		Nonlinear analysis time step

* **drift_capacity** - required for IDA, default to 10 (%)

		Assumed drift capacity for the building, beyond which the building is assumed to have collapsed

* **analysis_type** - required for All

		Analysis type to be run, list of strings, e.g. ["ST", "MA"] to run both ST and MA

* **system** - required for All, default to "space"

		May have two values:
			perimeter - exterior frames only as seismic lateral-load resisting frames
			space - all frames as seismic lateral-load resisting frames

* **hinge_model** - required for All, default to "Hysteretic"
		
		May have two values:
			Hysteretic - Hysteretic hinge models (uses offsets)
			Haselton - Haselton spring models (uses four node panel zones)

* **flag3d** - required for All, default to False

		False for 2D modelling
		True for 3D modelling

* **direction** - required for PO and ELF

		Direction of application for PO and ELF analysis
		0 stands for X direction, 1 stands for Y direction

* **export_at_each_step** - required for MSA and IDA, default to True

		True for exporting outputs for each time step (recommended)

* **period_assignment** - required for IDA, dictionary

		Period assignment ID for X and Y direction

* **periods_ida** - required for IDA

		List of float (periods) to use for IDA analysis

* **tcl_filename** - required for ST, MA, PO

		tcl filename necessary to generate tcl file models

	
</details>


### Examples
<details>
<a name="samples"></a>
<summary>Show/Hide</summary>
<br>

**3D building models**

Example 1: Static analysis - exampleStatic.py


Example 2: Modal analysis - exampleModel.py


Example 3: Static pushover analaysis - examplePushover.py


Example 4: MSA - exampleMSA.py


Example 5: IDA - exampleIDA.py


Example 6: Visualize - visualizeSPO.py


Example 7: 2D model with Haselton springs - exampleHaselton2D.py


</details>

### Future
<details>
<a name="future"></a>
<summary>Show/Hide</summary>
<br>

* [ ] Quality testing
* [ ] Haselton hinge models, example
* [x] 2D application examples
* [ ] Elastic models
* [x] 3D application examples
* [ ] Update solutionAlgorithm to incorporate interpolation functions (secondary analysis option, not recommended)


</details>
