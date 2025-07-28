from pymol import cmd

#restart
cmd.reinitialize()

# Show Seq window
#cmd.set("seq_view", 1)

# Fetch Papain
cmd.fetch ("1BP4")

# Hide everything
cmd.hide("all")

# Show everything in 30% grey and lines
cmd.show("lines", "all")
cmd.color("gray30", "all")

# Select cysteines and color them red
cmd.select("cys", "resn CYS")
cmd.color("red", "cys")

# Show disulfide bonds in cyan as sticks
cmd.select("ssbonds", "(resn CYS and name SG) within 2.2 of (resn CYS and name SG)")
cmd.show("sticks", "ssbonds")
cmd.color("cyan", "ssbonds")

# List of disulfide pairs from 1BP4 
ss_bonds = [
    ("/1BP4/A/A/CYS`22/SG", "/1BP4/A/A/CYS`63/SG"),
    ("/1BP4/A/A/CYS`56/SG", "/1BP4/A/A/CYS`95/SG"),
    ("/1BP4/A/A/CYS`153/SG", "/1BP4/A/A/CYS`200/SG")
]


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

# Label free cys: CYS 25 SG
cmd.label("/1BP4/A/A/CYS`25/SG", '"CYS25"')
cmd.show("spheres", "/1BP4/A/A/CYS`25/SG")
cmd.set("sphere_scale", 0.2, "/1BP4/A/A/CYS`25/SG")
cmd.color("green", "/1BP4/A/A/CYS`25/SG")


### Get distance ###
# Set up
from itertools import combinations
phantom_names = [f"phantom_SS_{i}" for i in range(1, 4)]
distance_threshold = 15


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
    #print(b)
    #print(sqrt(b / (8 * pi ** 2)))
    if b is None:
        return 0.0
    return sqrt(b / (8 * pi ** 2))

def get_phantom_uncertainty(atom1, atom2):
    sigma1 = get_atom_uncertainty(atom1)
    sigma2 = get_atom_uncertainty(atom2)
    return 0.5 * sqrt(sigma1**2 + sigma2**2)


# Loop through all combinations of pseudoatoms
for name1, name2 in combinations(phantom_names, 2):
    try:
        coord1 = cmd.get_atom_coords(name1)
        coord2 = cmd.get_atom_coords(name2)
        dist = ((sum((a - b) ** 2 for a, b in zip(coord1, coord2))) ** 0.5)
        
        if dist < distance_threshold:
            distance_label = f"dist_{name1}_{name2}"
            cmd.distance(distance_label, name1, name2)
            #cmd.set("dash_color", "cyan", distance_label)
            #cmd.set("dash_width", 2.0, distance_label)
    
            # Map back to original CYS atoms for uncertainty
            index1 = int(name1.split("_")[-1]) - 1
            index2 = int(name2.split("_")[-1]) - 1
            cys1_a, cys1_b = ss_bonds[index1]
            cys2_a, cys2_b = ss_bonds[index2]
            sigma1 = get_phantom_uncertainty(cys1_a, cys1_b)
            sigma2 = get_phantom_uncertainty(cys2_a, cys2_b)
            sigma_total = sqrt(sigma1**2 + sigma2**2)
            print(f"{distance_label}: {dist:.2f} ± {sigma_total:.2f} Å")
            #print ("hi", sigma1, sigma2)

    except:
        print(f"Could not compute or display distance between {name1} and {name2}")


# Compute distances from CYS25 to all phantom_SS_x and draw those under distance_threshold
cys25 = "/1BP4/A/A/CYS`25/SG"
for phantom in phantom_names:
    try:
        dist = cmd.get_distance(cys25, phantom)
        if dist < distance_threshold:
            label = f"dist_CYS25_{phantom}"
            cmd.distance(label, cys25, phantom)
            #cmd.set("dash_color", "magenta", label)
            #cmd.set("dash_width", 2.0, label)

            # Map back to original CYS atoms for uncertainty
            index1 = int(phantom.split("_")[-1]) - 1
            
            cys1_a, cys1_b = ss_bonds[index1]
            
            sigma1 = get_phantom_uncertainty(cys1_a, cys1_b)
            
            sigma2 = get_atom_uncertainty(cys25)
            
            sigma_total = sqrt(sigma1**2 + sigma2**2)
            print(f"{label}: {dist:.2f} ± {sigma_total:.2f} Å")
            #print("hii", sigma1, sigma2)
    except:
        print(f"Could not compute distance between CYS25 and {phantom}")
        
# Optional: zoom in
cmd.zoom("all")