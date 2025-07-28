from pymol import cmd

#restart
cmd.reinitialize()

# Show Seq window (disabled)
#cmd.set("seq_view", 1)

# Fetch BSA
cmd.fetch ("3v03")

# Set background color to white
cmd.bg_color("white")

# Hide everything
cmd.hide("all")

# Show everything in ribbon, color by chains, set ribbon transparency
cmd.show("ribbon","all")
util.color_chains("(3v03)",_self=cmd)
cmd.set("ribbon_transparency", 0.6)

# Select cysteines and color them red
cmd.select("cys", "resn CYS")

# Color Cys red (disabled)
#cmd.color("red", "cys")

# Show disulfide bonds in orange as sticks
cmd.select("ssbonds", "(resn CYS and name SG) within 2.2 of (resn CYS and name SG)")
cmd.show("sticks", "ssbonds")
cmd.color("orange", "ssbonds")
cmd.set("stick_radius", 0.7)

# List of disulfide pairs from 3v03
ss_bonds = [
    ("/3v03/A/A/CYS`53/SG", "/3v03/A/A/CYS`62/SG"),
    ("/3v03/A/A/CYS`75/SG", "/3v03/A/A/CYS`91/SG"),
    ("/3v03/A/A/CYS`90/SG", "/3v03/A/A/CYS`101/SG"),
    ("/3v03/A/A/CYS`123/SG", "/3v03/A/A/CYS`168/SG"),
    ("/3v03/A/A/CYS`167/SG", "/3v03/A/A/CYS`176/SG"),
    ("/3v03/A/A/CYS`199/SG", "/3v03/A/A/CYS`245/SG"),
    ("/3v03/A/A/CYS`244/SG", "/3v03/A/A/CYS`252/SG"),
    ("/3v03/A/A/CYS`264/SG", "/3v03/A/A/CYS`278/SG"),
    ("/3v03/A/A/CYS`277/SG", "/3v03/A/A/CYS`288/SG"),
    ("/3v03/A/A/CYS`315/SG", "/3v03/A/A/CYS`360/SG"),
    ("/3v03/A/A/CYS`359/SG", "/3v03/A/A/CYS`368/SG"),
    ("/3v03/A/A/CYS`391/SG", "/3v03/A/A/CYS`437/SG"),
    ("/3v03/A/A/CYS`436/SG", "/3v03/A/A/CYS`447/SG"),
    ("/3v03/A/A/CYS`460/SG", "/3v03/A/A/CYS`476/SG"),
    ("/3v03/A/A/CYS`475/SG", "/3v03/A/A/CYS`486/SG"),
    ("/3v03/A/A/CYS`513/SG", "/3v03/A/A/CYS`558/SG"),
    ("/3v03/A/A/CYS`557/SG", "/3v03/A/A/CYS`566/SG"),
    ("/3v03/B/B/CYS`53/SG", "/3v03/B/B/CYS`62/SG"),
    ("/3v03/B/B/CYS`75/SG", "/3v03/B/B/CYS`91/SG"),
    ("/3v03/B/B/CYS`90/SG", "/3v03/B/B/CYS`101/SG"),
    ("/3v03/B/B/CYS`123/SG", "/3v03/B/B/CYS`168/SG"),
    ("/3v03/B/B/CYS`167/SG", "/3v03/B/B/CYS`176/SG"),
    ("/3v03/B/B/CYS`199/SG", "/3v03/B/B/CYS`245/SG"),
    ("/3v03/B/B/CYS`244/SG", "/3v03/B/B/CYS`252/SG"),
    ("/3v03/B/B/CYS`264/SG", "/3v03/B/B/CYS`278/SG"),
    ("/3v03/B/B/CYS`277/SG", "/3v03/B/B/CYS`288/SG"),
    ("/3v03/B/B/CYS`315/SG", "/3v03/B/B/CYS`360/SG"),
    ("/3v03/B/B/CYS`359/SG", "/3v03/B/B/CYS`368/SG"),
    ("/3v03/B/B/CYS`391/SG", "/3v03/B/B/CYS`437/SG"),
    ("/3v03/B/B/CYS`436/SG", "/3v03/B/B/CYS`447/SG"),
    ("/3v03/B/B/CYS`460/SG", "/3v03/B/B/CYS`476/SG"),
    ("/3v03/B/B/CYS`475/SG", "/3v03/B/B/CYS`486/SG"),
    ("/3v03/B/B/CYS`513/SG", "/3v03/B/B/CYS`558/SG"),
    ("/3v03/B/B/CYS`557/SG", "/3v03/B/B/CYS`566/SG")
]

# Get Uncertainty_Basis
from math import sqrt, pi

def get_b_factor(selection):
    b_list = []
    cmd.iterate(selection, "b_list.append(b)", space={"b_list": b_list})
    if b_list:
        return b_list[0]  # assumes one atom in selection
    else:
        return None

def get_atom_uncertainty(atom_sel):
    b = get_b_factor(atom_sel)
    if b is None:
        return 0.0
    return sqrt(b / (8 * pi ** 2))

def get_phantom_uncertainty(atom1, atom2):
    sigma1 = get_atom_uncertainty(atom1)
    sigma2 = get_atom_uncertainty(atom2)
    return 0.5 * sqrt(sigma1**2 + sigma2**2)

# Make dummy atoms
for i, (cys1, cys2) in enumerate(ss_bonds, start=1):
    coord1 = cmd.get_atom_coords(cys1)
    coord2 = cmd.get_atom_coords(cys2)
    midpoint = [(a + b) / 2 for a, b in zip(coord1, coord2)]
    
    label = f"SS{i}"
    phantom_name = f"phantom_SS_{i}"
    cmd.pseudoatom(object=phantom_name, pos=midpoint, label=label, color="yellow")
    cmd.show("spheres", phantom_name)
    cmd.set("sphere_scale", 0.2, phantom_name)
    cmd.hide("labels", "all")
    

# Label free cys: CYS 34 SG in both chain
cmd.label("/3v03/A/A/CYS`34/SG", '"CYS_A_34"')
cmd.show("spheres", "/3v03/A/A/CYS`34/SG")
cmd.set("sphere_scale", 0.5, "/3v03/A/A/CYS`34/SG")
cmd.color("green", "/3v03/A/A/CYS`34/SG")

cmd.label("/3v03/B/B/CYS`34/SG", '"CYS_B_34"')
cmd.show("spheres", "/3v03/B/B/CYS`34/SG")
cmd.set("sphere_scale", 0.5, "/3v03/B/B/CYS`34/SG")
cmd.color("green", "/3v03/B/B/CYS`34/SG")


### Get distance ###
# Set up
from itertools import combinations
phantom_names = [f"phantom_SS_{i}" for i in range(1, 35)]

# Set distance threshold (Å) betweeen phantam atoms or free cys here:
distance_threshold = 15

# Loop through all combinations of phantom atoms
for name1, name2 in combinations(phantom_names, 2):
    try:
        coord1 = cmd.get_atom_coords(name1)
        coord2 = cmd.get_atom_coords(name2)
        dist = ((sum((a - b) ** 2 for a, b in zip(coord1, coord2))) ** 0.5)
        
        if dist < distance_threshold:
            distance_label = f"dist_{name1}_{name2}"
            cmd.distance(distance_label, name1, name2)
            cmd.set("dash_color", "magenta", distance_label)
            cmd.set("dash_width", 1.5, distance_label)

            # Map back to original CYS atoms for uncertainty
            index1 = int(name1.split("_")[-1]) - 1
            index2 = int(name2.split("_")[-1]) - 1
            cys1_a, cys1_b = ss_bonds[index1]
            cys2_a, cys2_b = ss_bonds[index2]
            sigma1 = get_phantom_uncertainty(cys1_a, cys1_b)
            sigma2 = get_phantom_uncertainty(cys2_a, cys2_b)
            sigma_total = sqrt(sigma1**2 + sigma2**2)
            print(f"{distance_label}: {dist:.2f} ± {sigma_total:.2f} Å")
    except:
        print(f"Could not compute or display distance between {name1} and {name2}")


### Compute distances from CYS34 to all phantom_SS_x and draw those under distance_threshold
# CYS34 in chain A
cys_A_34 = "/3v03/A/A/CYS`34/SG"
cmd.hide("labels", "/3v03/A/A/CYS`34/SG")
for phantom in phantom_names:
    try:
        dist = cmd.get_distance(cys_A_34, phantom)
        if dist < distance_threshold:
            label = f"dist_CYS_A_34_{phantom}"
            cmd.distance(label, cys_A_34, phantom)
            cmd.set("dash_color", "magenta", label)
            cmd.set("dash_width", 1.5, label)

            # Map back to original CYS atoms for uncertainty
            index1 = int(phantom.split("_")[-1]) - 1
            cys1_a, cys1_b = ss_bonds[index1]
            sigma1 = get_phantom_uncertainty(cys1_a, cys1_b)
            sigma2 = get_atom_uncertainty(cys_A_34)
            sigma_total = sqrt(sigma1**2 + sigma2**2)
            print(f"{label}: {dist:.2f} ± {sigma_total:.2f} Å")
            
    except:
        print(f"Could not compute distance between CYS_A_34 and {phantom}")

# CYS34 in chain B       
cys_B_34 = "/3v03/B/B/CYS`34/SG"
cmd.hide("labels", "/3v03/B/B/CYS`34/SG")
for phantom in phantom_names:
    try:
        dist = cmd.get_distance(cys_B_34, phantom)
        if dist < distance_threshold:
            label = f"dist_CYS_B_34_{phantom}"
            cmd.distance(label, cys_B_34, phantom)
            cmd.set("dash_color", "magenta", label)
            cmd.set("dash_width", 2.0, label)

            # Map back to original CYS atoms for uncertainty
            index1 = int(phantom.split("_")[-1]) - 1
            cys1_a, cys1_b = ss_bonds[index1]
            sigma1 = get_phantom_uncertainty(cys1_a, cys1_b)
            sigma2 = get_atom_uncertainty(cys_B_34)
            sigma_total = sqrt(sigma1**2 + sigma2**2)
            print(f"{label}: {dist:.2f} ± {sigma_total:.2f} Å")

    except:
        print(f"Could not compute distance between CYS_A_34 and {phantom}")
        
# label settings
cmd.hide("labels"    ,"All")
#cmd.label('''byca(cys)''', 'oneletter+resi')
#cmd.set("label_size", 14)
#cmd.set("label_shadow_mode", 1)

# View 1 (all)
view1 = (
        0.890873373,    0.300469279,   -0.340676486,\
        -0.374350935,    0.060814023,   -0.925290287,\
        -0.257304937,    0.951847136,    0.166658282,\
        0.000000000,    0.000000000, -465.119018555,\
        64.224395752,   25.798931122,   32.080112457,\
        334.588775635,  595.649291992,  -20.000000000 )

# View 2 (Zoom in_top, 4 ssbonds)
view2 = (
      0.617044866,    0.334960490,   -0.712076426,\
    -0.757793605,    0.008993667,   -0.652429163,\
    -0.212134540,    0.942184985,    0.259382099,\
     0.000036020,   -0.000593889,  -60.883811951,\
    71.414123535,   26.039836884,   54.976745605,\
    39.371078491,   82.400825500,  -20.000000000   )

# View 3 (Zoom in_bottom, cys34)
view3 = (
    0.945854843,    0.323722601,    0.023599001,\
     0.012369560,    0.036693495,   -0.999247909,\
    -0.324347585,    0.945437014,    0.030701196,\
    -0.000173680,   -0.000266965,  -75.748046875,\
    41.565170288,   15.416347504,   18.457775116,\
    64.188659668,   87.318565369,  -20.000000000)


# Choose view of interest here
cmd.set_view (view1)