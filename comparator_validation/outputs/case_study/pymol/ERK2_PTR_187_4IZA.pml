reinitialize
bg_color white
set cartoon_transparency, 0.65
set stick_radius, 0.16
set dash_radius, 0.06
set label_size, 16
set_color unmodified_gray80, [0.80, 0.80, 0.80]
set_color crystal_wheat_tint, [0.93, 0.82, 0.62]
set_color phosphofill_green_1, [0.00, 0.48, 0.18]
set_color phosphofill_green_2, [0.20, 0.65, 0.30]
set_color phosphofill_green_3, [0.50, 0.78, 0.45]
set_color pytms_blue_default, [0.25, 0.58, 0.95]
set_color pytms_blue_optimized, [0.02, 0.25, 0.75]
set_color ptmpsi_gold, [0.95, 0.62, 0.05]
set_color contact_crystal_green, [0.00, 0.48, 0.18]
set_color contact_maintained_gray, [0.45, 0.45, 0.45]
set_color contact_new_blue, [0.05, 0.35, 0.95]
set_color contact_lost_red, [0.88, 0.12, 0.10]
# ERK2 PTR 187 4IZA
# 5UMO author numbering places the unmodified ERK TEY Tyr at A:185; A:187 is Ala.
# Contact colors per prediction: maintained/recovered grey, new/extra blue, lost/missing red.
load "C:/Users/drpat/OneDrive - UGent/Desktop/PD/Phospho-fill/Benchmarkv2/Benchmark_combined_final/pdb_cache/5umo.cif", unmodified
load "C:/Users/drpat/OneDrive - UGent/Desktop/PD/Phospho-fill/Benchmarkv2/Benchmark_combined_final/pdb_cache/4iza.cif", phospho_crystal
align unmodified and chain A, phospho_crystal and chain A
hide everything, all
show cartoon, unmodified and chain A
show cartoon, phospho_crystal and chain A
hide everything, unmodified and not chain A
hide everything, phospho_crystal and not chain A
color unmodified_gray80, unmodified and chain A
color crystal_wheat_tint, phospho_crystal and chain A
select unmodified_site, unmodified and chain A and resi 185
select unmodified_marker, unmodified_site and name OH
select phospho_site, phospho_crystal and chain A and resi 187
select phospho_phosphate, phospho_site and name P+O1P+O2P+O3P+OP1+OP2+OP3
select phospho_phosphate_oxygen, phospho_site and name O1P+O2P+O3P+OP1+OP2+OP3
unbond (phospho_site and name P+O1P+O2P+O3P+OP1+OP2+OP3), (phospho_site)
bond (phospho_site and name OH), (phospho_site and name P)
bond (phospho_site and name P), (phospho_site and name O1P+OP1)
bond (phospho_site and name P), (phospho_site and name O2P+OP2)
bond (phospho_site and name P), (phospho_site and name O3P+OP3)
show sticks, unmodified_site or phospho_site
color crystal_wheat_tint, phospho_phosphate
select phospho_contacts_5A, byres (((phospho_crystal and chain A) within 5 of phospho_phosphate) and polymer.protein and not phospho_site)
show sticks, phospho_contacts_5A
select phospho_basic_residues_4A, byres (((phospho_crystal and chain A and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2) within 4 of phospho_phosphate_oxygen))
show sticks, phospho_basic_residues_4A
color contact_crystal_green, phospho_basic_residues_4A
distance phospho_basic_4A, phospho_phosphate_oxygen, ((phospho_crystal and chain A and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2) within 4 of phospho_phosphate_oxygen)
color contact_crystal_green, phospho_basic_4A
set dash_color, contact_crystal_green, phospho_basic_4A
set dash_width, 2.4, phospho_basic_4A
hide labels, phospho_basic_4A
select phospho_crystal_contact_residues, (phospho_crystal and chain A and resi 67+149+151+154+167+190)
load "C:/Users/drpat/OneDrive - UGent/Desktop/PD/Phospho-fill/Benchmarkv2/PhosphoFill-gihub/Benchmark/validation_structure_subset/outputs/archives/PTR/mode-joint__w-0.0__pack-0.2/P28482__4iza__A__target-187__ctx-187__mode-joint/stage3/phosphoFill1.cif", PhosphoFill1
align PhosphoFill1 and chain A, phospho_crystal and chain A
hide everything, PhosphoFill1
show cartoon, PhosphoFill1 and chain A
hide everything, PhosphoFill1 and not chain A
color phosphofill_green_1, PhosphoFill1 and chain A
select PhosphoFill1_site, PhosphoFill1 and chain A and resi 187
select PhosphoFill1_phosphate, PhosphoFill1_site and name P+O1P+O2P+O3P+OP1+OP2+OP3
unbond (PhosphoFill1_site and name P+O1P+O2P+O3P+OP1+OP2+OP3), (PhosphoFill1_site)
bond (PhosphoFill1_site and name OH), (PhosphoFill1_site and name P)
bond (PhosphoFill1_site and name P), (PhosphoFill1_site and name O1P+OP1)
bond (PhosphoFill1_site and name P), (PhosphoFill1_site and name O2P+OP2)
bond (PhosphoFill1_site and name P), (PhosphoFill1_site and name O3P+OP3)
show sticks, PhosphoFill1_site
select PhosphoFill1_contacts_5A, byres (((PhosphoFill1 and chain A) within 5 of PhosphoFill1_phosphate) and polymer.protein and not PhosphoFill1_site)
show sticks, PhosphoFill1_contacts_5A
select PhosphoFill1_maintained_basic_residues, byres ((PhosphoFill1 and chain A and resi 67+151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2))
show sticks, PhosphoFill1_maintained_basic_residues
color contact_maintained_gray, PhosphoFill1_maintained_basic_residues
distance PhosphoFill1_maintained_basic_contacts, PhosphoFill1_phosphate, ((PhosphoFill1 and chain A and resi 67+151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)), 4.0
color contact_maintained_gray, PhosphoFill1_maintained_basic_contacts
set dash_color, contact_maintained_gray, PhosphoFill1_maintained_basic_contacts
set dash_width, 2.2, PhosphoFill1_maintained_basic_contacts
hide labels, PhosphoFill1_maintained_basic_contacts
select PhosphoFill1_new_basic_residues, none
select PhosphoFill1_new_basic_contacts, none
select PhosphoFill1_lost_basic_residues, none
select PhosphoFill1_lost_basic_contacts, none
select PhosphoFill1_maintained_basic_atomic_donor_1, (PhosphoFill1 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1)
select PhosphoFill1_maintained_basic_atomic_phos_1, (PhosphoFill1_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PhosphoFill1_maintained_basic_atomic_contact_1, PhosphoFill1_maintained_basic_atomic_phos_1, PhosphoFill1_maintained_basic_atomic_donor_1, 4.0
color contact_maintained_gray, PhosphoFill1_maintained_basic_atomic_contact_1
set dash_color, contact_maintained_gray, PhosphoFill1_maintained_basic_atomic_contact_1
set dash_width, 2.0, PhosphoFill1_maintained_basic_atomic_contact_1
hide labels, PhosphoFill1_maintained_basic_atomic_contact_1
hide dashes, PhosphoFill1_maintained_basic_atomic_contact_1
select PhosphoFill1_maintained_basic_atomic_donor_2, (PhosphoFill1 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2)
select PhosphoFill1_maintained_basic_atomic_phos_2, (PhosphoFill1_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PhosphoFill1_maintained_basic_atomic_contact_2, PhosphoFill1_maintained_basic_atomic_phos_2, PhosphoFill1_maintained_basic_atomic_donor_2, 4.0
color contact_maintained_gray, PhosphoFill1_maintained_basic_atomic_contact_2
set dash_color, contact_maintained_gray, PhosphoFill1_maintained_basic_atomic_contact_2
set dash_width, 2.0, PhosphoFill1_maintained_basic_atomic_contact_2
hide labels, PhosphoFill1_maintained_basic_atomic_contact_2
hide dashes, PhosphoFill1_maintained_basic_atomic_contact_2
select PhosphoFill1_maintained_basic_atomic_donor_3, (PhosphoFill1 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)
select PhosphoFill1_maintained_basic_atomic_phos_3, (PhosphoFill1_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PhosphoFill1_maintained_basic_atomic_contact_3, PhosphoFill1_maintained_basic_atomic_phos_3, PhosphoFill1_maintained_basic_atomic_donor_3, 4.0
color contact_maintained_gray, PhosphoFill1_maintained_basic_atomic_contact_3
set dash_color, contact_maintained_gray, PhosphoFill1_maintained_basic_atomic_contact_3
set dash_width, 2.0, PhosphoFill1_maintained_basic_atomic_contact_3
hide labels, PhosphoFill1_maintained_basic_atomic_contact_3
hide dashes, PhosphoFill1_maintained_basic_atomic_contact_3
select PhosphoFill1_maintained_basic_atomic_residues, byres ((PhosphoFill1 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PhosphoFill1 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2) or (PhosphoFill1 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ))
select PhosphoFill1_maintained_basic_atomic_atoms, ((PhosphoFill1 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PhosphoFill1 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2) or (PhosphoFill1 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)) or (PhosphoFill1_phosphate)
group PhosphoFill1_maintained_basic_atomic_contacts, PhosphoFill1_maintained_basic_atomic_contact_1 PhosphoFill1_maintained_basic_atomic_contact_2 PhosphoFill1_maintained_basic_atomic_contact_3
select PhosphoFill1_new_basic_atomic_residues, none
select PhosphoFill1_new_basic_atomic_atoms, none
select PhosphoFill1_new_basic_atomic_contacts, none
select PhosphoFill1_lost_basic_atomic_residues, none
select PhosphoFill1_lost_basic_atomic_atoms, none
select PhosphoFill1_lost_basic_atomic_contacts, none
select PhosphoFill1_crystal_contact_residues, (PhosphoFill1 and chain A and resi 67+149+151+154+167+190)
load "C:/Users/drpat/OneDrive - UGent/Desktop/PD/Phospho-fill/Benchmarkv2/PhosphoFill-gihub/Benchmark/validation_structure_subset/outputs/archives/PTR/mode-joint__w-0.0__pack-0.2/P28482__4iza__A__target-187__ctx-187__mode-joint/stage3/phosphoFill2.cif", PhosphoFill2
align PhosphoFill2 and chain A, phospho_crystal and chain A
hide everything, PhosphoFill2
show cartoon, PhosphoFill2 and chain A
hide everything, PhosphoFill2 and not chain A
color phosphofill_green_2, PhosphoFill2 and chain A
select PhosphoFill2_site, PhosphoFill2 and chain A and resi 187
select PhosphoFill2_phosphate, PhosphoFill2_site and name P+O1P+O2P+O3P+OP1+OP2+OP3
unbond (PhosphoFill2_site and name P+O1P+O2P+O3P+OP1+OP2+OP3), (PhosphoFill2_site)
bond (PhosphoFill2_site and name OH), (PhosphoFill2_site and name P)
bond (PhosphoFill2_site and name P), (PhosphoFill2_site and name O1P+OP1)
bond (PhosphoFill2_site and name P), (PhosphoFill2_site and name O2P+OP2)
bond (PhosphoFill2_site and name P), (PhosphoFill2_site and name O3P+OP3)
show sticks, PhosphoFill2_site
select PhosphoFill2_contacts_5A, byres (((PhosphoFill2 and chain A) within 5 of PhosphoFill2_phosphate) and polymer.protein and not PhosphoFill2_site)
show sticks, PhosphoFill2_contacts_5A
select PhosphoFill2_maintained_basic_residues, byres ((PhosphoFill2 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2))
show sticks, PhosphoFill2_maintained_basic_residues
color contact_maintained_gray, PhosphoFill2_maintained_basic_residues
distance PhosphoFill2_maintained_basic_contacts, PhosphoFill2_phosphate, ((PhosphoFill2 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)), 4.0
color contact_maintained_gray, PhosphoFill2_maintained_basic_contacts
set dash_color, contact_maintained_gray, PhosphoFill2_maintained_basic_contacts
set dash_width, 2.2, PhosphoFill2_maintained_basic_contacts
hide labels, PhosphoFill2_maintained_basic_contacts
select PhosphoFill2_new_basic_residues, none
select PhosphoFill2_new_basic_contacts, none
select PhosphoFill2_lost_basic_residues, byres ((PhosphoFill2 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2))
show sticks, PhosphoFill2_lost_basic_residues
color contact_lost_red, PhosphoFill2_lost_basic_residues
distance PhosphoFill2_lost_basic_contacts, PhosphoFill2_phosphate, ((PhosphoFill2 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)), 12.0
color contact_lost_red, PhosphoFill2_lost_basic_contacts
set dash_color, contact_lost_red, PhosphoFill2_lost_basic_contacts
set dash_width, 2.2, PhosphoFill2_lost_basic_contacts
hide labels, PhosphoFill2_lost_basic_contacts
select PhosphoFill2_maintained_basic_atomic_donor_1, (PhosphoFill2 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1)
select PhosphoFill2_maintained_basic_atomic_phos_1, (PhosphoFill2_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PhosphoFill2_maintained_basic_atomic_contact_1, PhosphoFill2_maintained_basic_atomic_phos_1, PhosphoFill2_maintained_basic_atomic_donor_1, 4.0
color contact_maintained_gray, PhosphoFill2_maintained_basic_atomic_contact_1
set dash_color, contact_maintained_gray, PhosphoFill2_maintained_basic_atomic_contact_1
set dash_width, 2.0, PhosphoFill2_maintained_basic_atomic_contact_1
hide labels, PhosphoFill2_maintained_basic_atomic_contact_1
hide dashes, PhosphoFill2_maintained_basic_atomic_contact_1
select PhosphoFill2_maintained_basic_atomic_donor_2, (PhosphoFill2 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2)
select PhosphoFill2_maintained_basic_atomic_phos_2, (PhosphoFill2_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PhosphoFill2_maintained_basic_atomic_contact_2, PhosphoFill2_maintained_basic_atomic_phos_2, PhosphoFill2_maintained_basic_atomic_donor_2, 4.0
color contact_maintained_gray, PhosphoFill2_maintained_basic_atomic_contact_2
set dash_color, contact_maintained_gray, PhosphoFill2_maintained_basic_atomic_contact_2
set dash_width, 2.0, PhosphoFill2_maintained_basic_atomic_contact_2
hide labels, PhosphoFill2_maintained_basic_atomic_contact_2
hide dashes, PhosphoFill2_maintained_basic_atomic_contact_2
select PhosphoFill2_maintained_basic_atomic_residues, byres ((PhosphoFill2 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PhosphoFill2 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2))
select PhosphoFill2_maintained_basic_atomic_atoms, ((PhosphoFill2 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PhosphoFill2 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2)) or (PhosphoFill2_phosphate)
group PhosphoFill2_maintained_basic_atomic_contacts, PhosphoFill2_maintained_basic_atomic_contact_1 PhosphoFill2_maintained_basic_atomic_contact_2
select PhosphoFill2_new_basic_atomic_residues, none
select PhosphoFill2_new_basic_atomic_atoms, none
select PhosphoFill2_new_basic_atomic_contacts, none
select PhosphoFill2_lost_basic_atomic_donor_1, (PhosphoFill2 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)
select PhosphoFill2_lost_basic_atomic_phos_1, (PhosphoFill2_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PhosphoFill2_lost_basic_atomic_contact_1, PhosphoFill2_lost_basic_atomic_phos_1, PhosphoFill2_lost_basic_atomic_donor_1, 12.0
color contact_lost_red, PhosphoFill2_lost_basic_atomic_contact_1
set dash_color, contact_lost_red, PhosphoFill2_lost_basic_atomic_contact_1
set dash_width, 2.0, PhosphoFill2_lost_basic_atomic_contact_1
hide labels, PhosphoFill2_lost_basic_atomic_contact_1
hide dashes, PhosphoFill2_lost_basic_atomic_contact_1
select PhosphoFill2_lost_basic_atomic_residues, byres ((PhosphoFill2 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ))
select PhosphoFill2_lost_basic_atomic_atoms, ((PhosphoFill2 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)) or (PhosphoFill2_phosphate)
group PhosphoFill2_lost_basic_atomic_contacts, PhosphoFill2_lost_basic_atomic_contact_1
select PhosphoFill2_crystal_contact_residues, (PhosphoFill2 and chain A and resi 67+149+151+154+167+190)
load "C:/Users/drpat/OneDrive - UGent/Desktop/PD/Phospho-fill/Benchmarkv2/PhosphoFill-gihub/Benchmark/validation_structure_subset/outputs/archives/PTR/mode-joint__w-0.0__pack-0.2/P28482__4iza__A__target-187__ctx-187__mode-joint/stage3/phosphoFill3.cif", PhosphoFill3
align PhosphoFill3 and chain A, phospho_crystal and chain A
hide everything, PhosphoFill3
show cartoon, PhosphoFill3 and chain A
hide everything, PhosphoFill3 and not chain A
color phosphofill_green_3, PhosphoFill3 and chain A
select PhosphoFill3_site, PhosphoFill3 and chain A and resi 187
select PhosphoFill3_phosphate, PhosphoFill3_site and name P+O1P+O2P+O3P+OP1+OP2+OP3
unbond (PhosphoFill3_site and name P+O1P+O2P+O3P+OP1+OP2+OP3), (PhosphoFill3_site)
bond (PhosphoFill3_site and name OH), (PhosphoFill3_site and name P)
bond (PhosphoFill3_site and name P), (PhosphoFill3_site and name O1P+OP1)
bond (PhosphoFill3_site and name P), (PhosphoFill3_site and name O2P+OP2)
bond (PhosphoFill3_site and name P), (PhosphoFill3_site and name O3P+OP3)
show sticks, PhosphoFill3_site
select PhosphoFill3_contacts_5A, byres (((PhosphoFill3 and chain A) within 5 of PhosphoFill3_phosphate) and polymer.protein and not PhosphoFill3_site)
show sticks, PhosphoFill3_contacts_5A
select PhosphoFill3_maintained_basic_residues, byres ((PhosphoFill3 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2))
show sticks, PhosphoFill3_maintained_basic_residues
color contact_maintained_gray, PhosphoFill3_maintained_basic_residues
distance PhosphoFill3_maintained_basic_contacts, PhosphoFill3_phosphate, ((PhosphoFill3 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)), 4.0
color contact_maintained_gray, PhosphoFill3_maintained_basic_contacts
set dash_color, contact_maintained_gray, PhosphoFill3_maintained_basic_contacts
set dash_width, 2.2, PhosphoFill3_maintained_basic_contacts
hide labels, PhosphoFill3_maintained_basic_contacts
select PhosphoFill3_new_basic_residues, none
select PhosphoFill3_new_basic_contacts, none
select PhosphoFill3_lost_basic_residues, byres ((PhosphoFill3 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2))
show sticks, PhosphoFill3_lost_basic_residues
color contact_lost_red, PhosphoFill3_lost_basic_residues
distance PhosphoFill3_lost_basic_contacts, PhosphoFill3_phosphate, ((PhosphoFill3 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)), 12.0
color contact_lost_red, PhosphoFill3_lost_basic_contacts
set dash_color, contact_lost_red, PhosphoFill3_lost_basic_contacts
set dash_width, 2.2, PhosphoFill3_lost_basic_contacts
hide labels, PhosphoFill3_lost_basic_contacts
select PhosphoFill3_maintained_basic_atomic_donor_1, (PhosphoFill3 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1)
select PhosphoFill3_maintained_basic_atomic_phos_1, (PhosphoFill3_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PhosphoFill3_maintained_basic_atomic_contact_1, PhosphoFill3_maintained_basic_atomic_phos_1, PhosphoFill3_maintained_basic_atomic_donor_1, 4.0
color contact_maintained_gray, PhosphoFill3_maintained_basic_atomic_contact_1
set dash_color, contact_maintained_gray, PhosphoFill3_maintained_basic_atomic_contact_1
set dash_width, 2.0, PhosphoFill3_maintained_basic_atomic_contact_1
hide labels, PhosphoFill3_maintained_basic_atomic_contact_1
hide dashes, PhosphoFill3_maintained_basic_atomic_contact_1
select PhosphoFill3_maintained_basic_atomic_donor_2, (PhosphoFill3 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2)
select PhosphoFill3_maintained_basic_atomic_phos_2, (PhosphoFill3_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PhosphoFill3_maintained_basic_atomic_contact_2, PhosphoFill3_maintained_basic_atomic_phos_2, PhosphoFill3_maintained_basic_atomic_donor_2, 4.0
color contact_maintained_gray, PhosphoFill3_maintained_basic_atomic_contact_2
set dash_color, contact_maintained_gray, PhosphoFill3_maintained_basic_atomic_contact_2
set dash_width, 2.0, PhosphoFill3_maintained_basic_atomic_contact_2
hide labels, PhosphoFill3_maintained_basic_atomic_contact_2
hide dashes, PhosphoFill3_maintained_basic_atomic_contact_2
select PhosphoFill3_maintained_basic_atomic_residues, byres ((PhosphoFill3 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PhosphoFill3 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2))
select PhosphoFill3_maintained_basic_atomic_atoms, ((PhosphoFill3 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PhosphoFill3 and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2)) or (PhosphoFill3_phosphate)
group PhosphoFill3_maintained_basic_atomic_contacts, PhosphoFill3_maintained_basic_atomic_contact_1 PhosphoFill3_maintained_basic_atomic_contact_2
select PhosphoFill3_new_basic_atomic_residues, none
select PhosphoFill3_new_basic_atomic_atoms, none
select PhosphoFill3_new_basic_atomic_contacts, none
select PhosphoFill3_lost_basic_atomic_donor_1, (PhosphoFill3 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)
select PhosphoFill3_lost_basic_atomic_phos_1, (PhosphoFill3_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PhosphoFill3_lost_basic_atomic_contact_1, PhosphoFill3_lost_basic_atomic_phos_1, PhosphoFill3_lost_basic_atomic_donor_1, 12.0
color contact_lost_red, PhosphoFill3_lost_basic_atomic_contact_1
set dash_color, contact_lost_red, PhosphoFill3_lost_basic_atomic_contact_1
set dash_width, 2.0, PhosphoFill3_lost_basic_atomic_contact_1
hide labels, PhosphoFill3_lost_basic_atomic_contact_1
hide dashes, PhosphoFill3_lost_basic_atomic_contact_1
select PhosphoFill3_lost_basic_atomic_residues, byres ((PhosphoFill3 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ))
select PhosphoFill3_lost_basic_atomic_atoms, ((PhosphoFill3 and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)) or (PhosphoFill3_phosphate)
group PhosphoFill3_lost_basic_atomic_contacts, PhosphoFill3_lost_basic_atomic_contact_1
select PhosphoFill3_crystal_contact_residues, (PhosphoFill3 and chain A and resi 67+149+151+154+167+190)
load "C:/Users/drpat/OneDrive - UGent/Desktop/PD/Phospho-fill/Benchmarkv2/Benchmark_combined_final/case_study_outputs/pytms_work_ERK2_Y187/P28482_PTR_4iza_A_187_pytms.pdb", PyTMs_default
align PyTMs_default and chain A, phospho_crystal and chain A
hide everything, PyTMs_default
show cartoon, PyTMs_default and chain A
hide everything, PyTMs_default and not chain A
color pytms_blue_default, PyTMs_default and chain A
select PyTMs_default_site, PyTMs_default and chain A and resi 187
select PyTMs_default_phosphate, PyTMs_default_site and name P+O1P+O2P+O3P+OP1+OP2+OP3
unbond (PyTMs_default_site and name P+O1P+O2P+O3P+OP1+OP2+OP3), (PyTMs_default_site)
bond (PyTMs_default_site and name OH), (PyTMs_default_site and name P)
bond (PyTMs_default_site and name P), (PyTMs_default_site and name O1P+OP1)
bond (PyTMs_default_site and name P), (PyTMs_default_site and name O2P+OP2)
bond (PyTMs_default_site and name P), (PyTMs_default_site and name O3P+OP3)
show sticks, PyTMs_default_site
select PyTMs_default_contacts_5A, byres (((PyTMs_default and chain A) within 5 of PyTMs_default_phosphate) and polymer.protein and not PyTMs_default_site)
show sticks, PyTMs_default_contacts_5A
select PyTMs_default_maintained_basic_residues, byres ((PyTMs_default and chain A and resi 67+151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2))
show sticks, PyTMs_default_maintained_basic_residues
color contact_maintained_gray, PyTMs_default_maintained_basic_residues
distance PyTMs_default_maintained_basic_contacts, PyTMs_default_phosphate, ((PyTMs_default and chain A and resi 67+151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)), 4.0
color contact_maintained_gray, PyTMs_default_maintained_basic_contacts
set dash_color, contact_maintained_gray, PyTMs_default_maintained_basic_contacts
set dash_width, 2.2, PyTMs_default_maintained_basic_contacts
hide labels, PyTMs_default_maintained_basic_contacts
select PyTMs_default_new_basic_residues, none
select PyTMs_default_new_basic_contacts, none
select PyTMs_default_lost_basic_residues, none
select PyTMs_default_lost_basic_contacts, none
select PyTMs_default_maintained_basic_atomic_donor_1, (PyTMs_default and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1)
select PyTMs_default_maintained_basic_atomic_phos_1, (PyTMs_default_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PyTMs_default_maintained_basic_atomic_contact_1, PyTMs_default_maintained_basic_atomic_phos_1, PyTMs_default_maintained_basic_atomic_donor_1, 4.0
color contact_maintained_gray, PyTMs_default_maintained_basic_atomic_contact_1
set dash_color, contact_maintained_gray, PyTMs_default_maintained_basic_atomic_contact_1
set dash_width, 2.0, PyTMs_default_maintained_basic_atomic_contact_1
hide labels, PyTMs_default_maintained_basic_atomic_contact_1
hide dashes, PyTMs_default_maintained_basic_atomic_contact_1
select PyTMs_default_maintained_basic_atomic_donor_2, (PyTMs_default and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2)
select PyTMs_default_maintained_basic_atomic_phos_2, (PyTMs_default_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PyTMs_default_maintained_basic_atomic_contact_2, PyTMs_default_maintained_basic_atomic_phos_2, PyTMs_default_maintained_basic_atomic_donor_2, 4.0
color contact_maintained_gray, PyTMs_default_maintained_basic_atomic_contact_2
set dash_color, contact_maintained_gray, PyTMs_default_maintained_basic_atomic_contact_2
set dash_width, 2.0, PyTMs_default_maintained_basic_atomic_contact_2
hide labels, PyTMs_default_maintained_basic_atomic_contact_2
hide dashes, PyTMs_default_maintained_basic_atomic_contact_2
select PyTMs_default_maintained_basic_atomic_donor_3, (PyTMs_default and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)
select PyTMs_default_maintained_basic_atomic_phos_3, (PyTMs_default_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PyTMs_default_maintained_basic_atomic_contact_3, PyTMs_default_maintained_basic_atomic_phos_3, PyTMs_default_maintained_basic_atomic_donor_3, 4.0
color contact_maintained_gray, PyTMs_default_maintained_basic_atomic_contact_3
set dash_color, contact_maintained_gray, PyTMs_default_maintained_basic_atomic_contact_3
set dash_width, 2.0, PyTMs_default_maintained_basic_atomic_contact_3
hide labels, PyTMs_default_maintained_basic_atomic_contact_3
hide dashes, PyTMs_default_maintained_basic_atomic_contact_3
select PyTMs_default_maintained_basic_atomic_residues, byres ((PyTMs_default and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PyTMs_default and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2) or (PyTMs_default and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ))
select PyTMs_default_maintained_basic_atomic_atoms, ((PyTMs_default and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PyTMs_default and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2) or (PyTMs_default and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)) or (PyTMs_default_phosphate)
group PyTMs_default_maintained_basic_atomic_contacts, PyTMs_default_maintained_basic_atomic_contact_1 PyTMs_default_maintained_basic_atomic_contact_2 PyTMs_default_maintained_basic_atomic_contact_3
select PyTMs_default_new_basic_atomic_residues, none
select PyTMs_default_new_basic_atomic_atoms, none
select PyTMs_default_new_basic_atomic_contacts, none
select PyTMs_default_lost_basic_atomic_residues, none
select PyTMs_default_lost_basic_atomic_atoms, none
select PyTMs_default_lost_basic_atomic_contacts, none
select PyTMs_default_crystal_contact_residues, (PyTMs_default and chain A and resi 67+149+151+154+167+190)
load "C:/Users/drpat/OneDrive - UGent/Desktop/PD/Phospho-fill/Benchmarkv2/Benchmark_combined_final/case_study_outputs/pytms_rotation_work_ERK2_Y187/P28482_PTR_4iza_A_187_pytms.pdb", PyTMs_optimized
align PyTMs_optimized and chain A, phospho_crystal and chain A
hide everything, PyTMs_optimized
show cartoon, PyTMs_optimized and chain A
hide everything, PyTMs_optimized and not chain A
color pytms_blue_optimized, PyTMs_optimized and chain A
select PyTMs_optimized_site, PyTMs_optimized and chain A and resi 187
select PyTMs_optimized_phosphate, PyTMs_optimized_site and name P+O1P+O2P+O3P+OP1+OP2+OP3
unbond (PyTMs_optimized_site and name P+O1P+O2P+O3P+OP1+OP2+OP3), (PyTMs_optimized_site)
bond (PyTMs_optimized_site and name OH), (PyTMs_optimized_site and name P)
bond (PyTMs_optimized_site and name P), (PyTMs_optimized_site and name O1P+OP1)
bond (PyTMs_optimized_site and name P), (PyTMs_optimized_site and name O2P+OP2)
bond (PyTMs_optimized_site and name P), (PyTMs_optimized_site and name O3P+OP3)
show sticks, PyTMs_optimized_site
select PyTMs_optimized_contacts_5A, byres (((PyTMs_optimized and chain A) within 5 of PyTMs_optimized_phosphate) and polymer.protein and not PyTMs_optimized_site)
show sticks, PyTMs_optimized_contacts_5A
select PyTMs_optimized_maintained_basic_residues, byres ((PyTMs_optimized and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2))
show sticks, PyTMs_optimized_maintained_basic_residues
color contact_maintained_gray, PyTMs_optimized_maintained_basic_residues
distance PyTMs_optimized_maintained_basic_contacts, PyTMs_optimized_phosphate, ((PyTMs_optimized and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)), 4.0
color contact_maintained_gray, PyTMs_optimized_maintained_basic_contacts
set dash_color, contact_maintained_gray, PyTMs_optimized_maintained_basic_contacts
set dash_width, 2.2, PyTMs_optimized_maintained_basic_contacts
hide labels, PyTMs_optimized_maintained_basic_contacts
select PyTMs_optimized_new_basic_residues, none
select PyTMs_optimized_new_basic_contacts, none
select PyTMs_optimized_lost_basic_residues, byres ((PyTMs_optimized and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2))
show sticks, PyTMs_optimized_lost_basic_residues
color contact_lost_red, PyTMs_optimized_lost_basic_residues
distance PyTMs_optimized_lost_basic_contacts, PyTMs_optimized_phosphate, ((PyTMs_optimized and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)), 12.0
color contact_lost_red, PyTMs_optimized_lost_basic_contacts
set dash_color, contact_lost_red, PyTMs_optimized_lost_basic_contacts
set dash_width, 2.2, PyTMs_optimized_lost_basic_contacts
hide labels, PyTMs_optimized_lost_basic_contacts
select PyTMs_optimized_maintained_basic_atomic_donor_1, (PyTMs_optimized and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)
select PyTMs_optimized_maintained_basic_atomic_phos_1, (PyTMs_optimized_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PyTMs_optimized_maintained_basic_atomic_contact_1, PyTMs_optimized_maintained_basic_atomic_phos_1, PyTMs_optimized_maintained_basic_atomic_donor_1, 4.0
color contact_maintained_gray, PyTMs_optimized_maintained_basic_atomic_contact_1
set dash_color, contact_maintained_gray, PyTMs_optimized_maintained_basic_atomic_contact_1
set dash_width, 2.0, PyTMs_optimized_maintained_basic_atomic_contact_1
hide labels, PyTMs_optimized_maintained_basic_atomic_contact_1
hide dashes, PyTMs_optimized_maintained_basic_atomic_contact_1
select PyTMs_optimized_maintained_basic_atomic_residues, byres ((PyTMs_optimized and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ))
select PyTMs_optimized_maintained_basic_atomic_atoms, ((PyTMs_optimized and chain A and resi 151 and polymer.protein and resn ARG+LYS+HIS and name NZ)) or (PyTMs_optimized_phosphate)
group PyTMs_optimized_maintained_basic_atomic_contacts, PyTMs_optimized_maintained_basic_atomic_contact_1
select PyTMs_optimized_new_basic_atomic_residues, none
select PyTMs_optimized_new_basic_atomic_atoms, none
select PyTMs_optimized_new_basic_atomic_contacts, none
select PyTMs_optimized_lost_basic_atomic_donor_1, (PyTMs_optimized and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1)
select PyTMs_optimized_lost_basic_atomic_phos_1, (PyTMs_optimized_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PyTMs_optimized_lost_basic_atomic_contact_1, PyTMs_optimized_lost_basic_atomic_phos_1, PyTMs_optimized_lost_basic_atomic_donor_1, 12.0
color contact_lost_red, PyTMs_optimized_lost_basic_atomic_contact_1
set dash_color, contact_lost_red, PyTMs_optimized_lost_basic_atomic_contact_1
set dash_width, 2.0, PyTMs_optimized_lost_basic_atomic_contact_1
hide labels, PyTMs_optimized_lost_basic_atomic_contact_1
hide dashes, PyTMs_optimized_lost_basic_atomic_contact_1
select PyTMs_optimized_lost_basic_atomic_donor_2, (PyTMs_optimized and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2)
select PyTMs_optimized_lost_basic_atomic_phos_2, (PyTMs_optimized_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PyTMs_optimized_lost_basic_atomic_contact_2, PyTMs_optimized_lost_basic_atomic_phos_2, PyTMs_optimized_lost_basic_atomic_donor_2, 12.0
color contact_lost_red, PyTMs_optimized_lost_basic_atomic_contact_2
set dash_color, contact_lost_red, PyTMs_optimized_lost_basic_atomic_contact_2
set dash_width, 2.0, PyTMs_optimized_lost_basic_atomic_contact_2
hide labels, PyTMs_optimized_lost_basic_atomic_contact_2
hide dashes, PyTMs_optimized_lost_basic_atomic_contact_2
select PyTMs_optimized_lost_basic_atomic_residues, byres ((PyTMs_optimized and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PyTMs_optimized and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2))
select PyTMs_optimized_lost_basic_atomic_atoms, ((PyTMs_optimized and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PyTMs_optimized and chain A and resi 67 and polymer.protein and resn ARG+LYS+HIS and name NH2)) or (PyTMs_optimized_phosphate)
group PyTMs_optimized_lost_basic_atomic_contacts, PyTMs_optimized_lost_basic_atomic_contact_1 PyTMs_optimized_lost_basic_atomic_contact_2
select PyTMs_optimized_crystal_contact_residues, (PyTMs_optimized and chain A and resi 67+149+151+154+167+190)
load "C:/Users/drpat/OneDrive - UGent/Desktop/PD/Phospho-fill/Benchmarkv2/Benchmark_combined_final/case_study_outputs/ptmpsi_work_ERK2_Y187/P28482_PTR_4iza_A_187_ptmpsi.pdb", PTM_Psi
align PTM_Psi and chain A, phospho_crystal and chain A
hide everything, PTM_Psi
show cartoon, PTM_Psi and chain A
hide everything, PTM_Psi and not chain A
color ptmpsi_gold, PTM_Psi and chain A
select PTM_Psi_site, PTM_Psi and chain A and resi 176
select PTM_Psi_phosphate, PTM_Psi_site and name P+O1P+O2P+O3P+OP1+OP2+OP3
unbond (PTM_Psi_site and name P+O1P+O2P+O3P+OP1+OP2+OP3), (PTM_Psi_site)
bond (PTM_Psi_site and name OH), (PTM_Psi_site and name P)
bond (PTM_Psi_site and name P), (PTM_Psi_site and name O1P+OP1)
bond (PTM_Psi_site and name P), (PTM_Psi_site and name O2P+OP2)
bond (PTM_Psi_site and name P), (PTM_Psi_site and name O3P+OP3)
show sticks, PTM_Psi_site
select PTM_Psi_contacts_5A, byres (((PTM_Psi and chain A) within 5 of PTM_Psi_phosphate) and polymer.protein and not PTM_Psi_site)
show sticks, PTM_Psi_contacts_5A
select PTM_Psi_maintained_basic_residues, none
select PTM_Psi_maintained_basic_contacts, none
select PTM_Psi_new_basic_residues, none
select PTM_Psi_new_basic_contacts, none
select PTM_Psi_lost_basic_residues, byres ((PTM_Psi and chain A and resi 57+141 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2))
show sticks, PTM_Psi_lost_basic_residues
color contact_lost_red, PTM_Psi_lost_basic_residues
distance PTM_Psi_lost_basic_contacts, PTM_Psi_phosphate, ((PTM_Psi and chain A and resi 57+141 and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)), 12.0
color contact_lost_red, PTM_Psi_lost_basic_contacts
set dash_color, contact_lost_red, PTM_Psi_lost_basic_contacts
set dash_width, 2.2, PTM_Psi_lost_basic_contacts
hide labels, PTM_Psi_lost_basic_contacts
select PTM_Psi_maintained_basic_atomic_residues, none
select PTM_Psi_maintained_basic_atomic_atoms, none
select PTM_Psi_maintained_basic_atomic_contacts, none
select PTM_Psi_new_basic_atomic_residues, none
select PTM_Psi_new_basic_atomic_atoms, none
select PTM_Psi_new_basic_atomic_contacts, none
select PTM_Psi_lost_basic_atomic_donor_1, (PTM_Psi and chain A and resi 57 and polymer.protein and resn ARG+LYS+HIS and name NH1)
select PTM_Psi_lost_basic_atomic_phos_1, (PTM_Psi_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PTM_Psi_lost_basic_atomic_contact_1, PTM_Psi_lost_basic_atomic_phos_1, PTM_Psi_lost_basic_atomic_donor_1, 12.0
color contact_lost_red, PTM_Psi_lost_basic_atomic_contact_1
set dash_color, contact_lost_red, PTM_Psi_lost_basic_atomic_contact_1
set dash_width, 2.0, PTM_Psi_lost_basic_atomic_contact_1
hide labels, PTM_Psi_lost_basic_atomic_contact_1
hide dashes, PTM_Psi_lost_basic_atomic_contact_1
select PTM_Psi_lost_basic_atomic_donor_2, (PTM_Psi and chain A and resi 57 and polymer.protein and resn ARG+LYS+HIS and name NH2)
select PTM_Psi_lost_basic_atomic_phos_2, (PTM_Psi_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PTM_Psi_lost_basic_atomic_contact_2, PTM_Psi_lost_basic_atomic_phos_2, PTM_Psi_lost_basic_atomic_donor_2, 12.0
color contact_lost_red, PTM_Psi_lost_basic_atomic_contact_2
set dash_color, contact_lost_red, PTM_Psi_lost_basic_atomic_contact_2
set dash_width, 2.0, PTM_Psi_lost_basic_atomic_contact_2
hide labels, PTM_Psi_lost_basic_atomic_contact_2
hide dashes, PTM_Psi_lost_basic_atomic_contact_2
select PTM_Psi_lost_basic_atomic_donor_3, (PTM_Psi and chain A and resi 141 and polymer.protein and resn ARG+LYS+HIS and name NZ)
select PTM_Psi_lost_basic_atomic_phos_3, (PTM_Psi_phosphate and name O1P+O2P+O3P+OP1+OP2+OP3)
distance PTM_Psi_lost_basic_atomic_contact_3, PTM_Psi_lost_basic_atomic_phos_3, PTM_Psi_lost_basic_atomic_donor_3, 12.0
color contact_lost_red, PTM_Psi_lost_basic_atomic_contact_3
set dash_color, contact_lost_red, PTM_Psi_lost_basic_atomic_contact_3
set dash_width, 2.0, PTM_Psi_lost_basic_atomic_contact_3
hide labels, PTM_Psi_lost_basic_atomic_contact_3
hide dashes, PTM_Psi_lost_basic_atomic_contact_3
select PTM_Psi_lost_basic_atomic_residues, byres ((PTM_Psi and chain A and resi 57 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PTM_Psi and chain A and resi 57 and polymer.protein and resn ARG+LYS+HIS and name NH2) or (PTM_Psi and chain A and resi 141 and polymer.protein and resn ARG+LYS+HIS and name NZ))
select PTM_Psi_lost_basic_atomic_atoms, ((PTM_Psi and chain A and resi 57 and polymer.protein and resn ARG+LYS+HIS and name NH1) or (PTM_Psi and chain A and resi 57 and polymer.protein and resn ARG+LYS+HIS and name NH2) or (PTM_Psi and chain A and resi 141 and polymer.protein and resn ARG+LYS+HIS and name NZ)) or (PTM_Psi_phosphate)
group PTM_Psi_lost_basic_atomic_contacts, PTM_Psi_lost_basic_atomic_contact_1 PTM_Psi_lost_basic_atomic_contact_2 PTM_Psi_lost_basic_atomic_contact_3
select PTM_Psi_crystal_contact_residues, (PTM_Psi and chain A and resi 57+139+141+144+157+178)
# Convenience groups: show dashes, PyTMs_maintained_basic_contacts / PyTMs_new_basic_contacts / PyTMs_lost_basic_contacts
# Donor-atom groups: show dashes, PyTMs_maintained_basic_atomic_contacts / PyTMs_new_basic_atomic_contacts / PyTMs_lost_basic_atomic_contacts
group PhosphoFill_maintained_basic_contacts, PhosphoFill1_maintained_basic_contacts PhosphoFill2_maintained_basic_contacts PhosphoFill3_maintained_basic_contacts
select PhosphoFill_new_basic_contacts, none
group PhosphoFill_lost_basic_contacts, PhosphoFill2_lost_basic_contacts PhosphoFill3_lost_basic_contacts
group PyTMs_maintained_basic_contacts, PyTMs_default_maintained_basic_contacts PyTMs_optimized_maintained_basic_contacts
select PyTMs_new_basic_contacts, none
group PyTMs_lost_basic_contacts, PyTMs_optimized_lost_basic_contacts
group PhosphoFill_maintained_basic_atomic_contacts, PhosphoFill1_maintained_basic_atomic_contacts PhosphoFill2_maintained_basic_atomic_contacts PhosphoFill3_maintained_basic_atomic_contacts
select PhosphoFill_new_basic_atomic_contacts, none
group PhosphoFill_lost_basic_atomic_contacts, PhosphoFill2_lost_basic_atomic_contacts PhosphoFill3_lost_basic_atomic_contacts
group PyTMs_maintained_basic_atomic_contacts, PyTMs_default_maintained_basic_atomic_contacts PyTMs_optimized_maintained_basic_atomic_contacts
select PyTMs_new_basic_atomic_contacts, none
group PyTMs_lost_basic_atomic_contacts, PyTMs_optimized_lost_basic_atomic_contacts
orient phospho_site
zoom phospho_site or phospho_contacts_5A, 12
set_view (\
    0.7, 0.1, 0.7,\
    0.0, 1.0, -0.2,\
    -0.7, 0.2, 0.7,\
    0.0, 0.0, -120.0,\
    0.0, 0.0, 0.0,\
    80.0, 160.0, -20.0 )
