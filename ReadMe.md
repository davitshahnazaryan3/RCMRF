<h1 align="center">Generic RC MRF model creator</h1>

Procedures to create a reinforced concrete building (3D), or a frame (2D) nonlinear structural model
consisting of moment resisting frames (MRF) as the primary lateral load resisting system.

**Required libraries**: requirements.txt


### Table of Contents
<details open>
<summary>Show/Hide</summary>
<br>

1. [Step-by-step procedure](#steps)

</details>

### Step-by-step procedure
<details>
<a name="steps"></a>
<summary>Show/Hide</summary>
<br>

* Master file
* Model file (Model)
	1. Define model dimensions and DOFs - create_model
	2. Create model nodes and defines boudnary conditions - create_nodes, client.Geometry
	3. Define transformations - define_transformations
	4. Define joint materials - joint_materials
	5. Define rotational springs - rot_springs, client.Sections
	6. Define bilinear springs - bilin_springs, client.Sections
	7. Create elements - create_elements
	8. Create joints - create_joints
	9. Define loads - define_loads
	10. Define P-Delta columns - define_pdelta_columns
	11. Define masses - define_masses
	12. Perform analysis - perform_analysis
		
			* ELF - equivalent lateral force - analysis.Static
			* ST - static gravity - analysis.Static
			* MA - modal - analysis.Modal
			* PO - static pushover - analysis.SPO
			* TH - time history/IDA - analysis.IDA_HTF, analysis.SolutionAlgorithm
		
	13. Set recorders - set_recorders, client.Recorders
	
</details>

