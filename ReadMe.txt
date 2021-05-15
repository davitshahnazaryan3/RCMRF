Generic RC MRF model creator

Includes both 3D and 2D modelling

Requires: OpenSeesPy

Step-by-step procedure
*A. Master file
B. Model file (Model)
	1. Defines model dimensions and DOFs - create_model
	2. Creates model nodes and defines boudnary conditions - create_nodes, client.Geometry
	3. Defines transformations - define_transformations
	4. Defines joint materials - joint_materials
	5. Defines rotational springs - rot_springs, client.Sections
	6. Defines bilinear springs - bilin_springs, client.Sections
	7. Creates elements - create_elements
	8. Creates joints - create_joints
	9. Defines loads - define_loads
	10. Defines P-Delta columns - define_pdelta_columns
	11. Defines masses - define_masses
	12. Performs analysis - perform_analysis
		a. ELF - equivalent lateral force - analysis.Static
		b. ST - static gravity - analysis.Static
		c. MA - modal - analysis.Modal
		d. PO - static pushover - analysis.SPO
		*e. TH - time history/IDA - analysis.IDA_HTF, analysis.SolutionAlgorithm
	13. Set recorders - set_recorders, client.Recorders



------ Issues to address
* Example record 29: for large IM a very small displacement is recorded. Find out why.


Future updates

- Space frames
- Haselton spring models for a 3D model











