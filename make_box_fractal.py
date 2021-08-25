import numpy as np
import sys
import maschberger_IMF

def CoM(n, x, m):

    #Find the centre of mass and move it to 0,0,0

    #Initialise
    xCoM, yCoM, zCoM, total_mass = 0,0,0,0
 
    #Go through each star
    for i in range(n):

        #Add the masses to the total
        total_mass = total_mass + m[i]

        #Add the weighted x, y, and z coordinates to their total
        xCoM = xCoM + ( x[0][i] * m[i] )
        yCoM = yCoM + ( x[1][i] * m[i] )
        zCoM = zCoM + ( x[2][i] * m[i] )


    #Find the CoM
    CoM = [ xCoM / total_mass, yCoM / total_mass, zCoM / total_mass  ]

        
    #Move all the coordinates, so that the CoM is centred on 0.
    for i in range(n):

        #Go through x, y, and z
        for j in range(3):

            #Shift the coordinated by the CoM amount
            x[j][i] = x[j][i] - CoM[j]


    return x;


def get_pos_and_vel(box_positions, box_vels, box, new_box_positions, new_box_vels, prob_survive, side_length, level, dx, dy, dz):

    # Given a child position decide if there's a star there and
    # give it a velocity

    # Get the coordinates within the box to look at creating a new star
    x_coord = box_positions[0][box] + (dx*side_length/2.)
    y_coord = box_positions[1][box] + (dy*side_length/2.)
    z_coord = box_positions[2][box] + (dz*side_length/2.)

    
    #Each one has an prob_survive chance of having a star in it
    #Chose a random number to decide
    if np.random.uniform() < prob_survive:

        success = True

        #There is a star, so give it a position and velocity
        #Add a bit of noise so it's no too griddy.
        noise_factor = 0.1*side_length*(1/3.)#abs(dx)
        new_box_positions[0].append(np.random.normal(scale=noise_factor, loc=x_coord))
        new_box_positions[1].append(np.random.normal(scale=noise_factor, loc=y_coord))
        new_box_positions[2].append(np.random.normal(scale=noise_factor, loc=z_coord))

        noise_factor = (1./level)**(2./3.)
        new_box_vels[0].append(np.random.normal(scale=noise_factor, loc=box_vels[0][box]))
        new_box_vels[1].append(np.random.normal(scale=noise_factor, loc=box_vels[1][box]))
        new_box_vels[2].append(np.random.normal(scale=noise_factor, loc=box_vels[2][box]))

    else:
        success = False

    return success, new_box_positions, new_box_vels

# This function makes a fractal where each child has a prob_survive
# probability of survival. If vel_structe is set to false then shuffle
# the stars velocities to remove any velocity structure before outputting.
def make_fractal(n_stars, prob_survive, vel_struct=True):

    #Setup the starting point
    box_positions = [[0.], [0.], [0.]]
    box_vels = [[0.], [0.], [0.]]
    side_length = 2.
    level = 0.

    #Create at least five times the desired number of stars
    while len(box_positions[0]) < 10*n_stars:
    
        level += 1
    
        new_box_positions = [[], [], []]
        new_box_vels = [[], [], []]

        for box in range(len(box_positions[0])):

            # Split the box into 27
            # Calculate the shift of each box from the centre
            offset = [-2/3., 0, 2/3.] 
            fail_list = []
            for dx in offset:
                for dy in offset:
                    for dz in offset:
                    
                        success, new_box_positions, new_box_vels = get_pos_and_vel(box_positions, box_vels, box, new_box_positions, new_box_vels, prob_survive, side_length, level, dx, dy, dz)

                        # Log the boxes where stars did not get created
                        if not success:
                            fail_list.append([dx, dy, dz])


            # If at a high level there happens to be a lot of infant moratlity just due
            # to stochasticity it'd lead to a lot more substructure than there should be
            # So if a certain number or more failed re-do them until there's enough.
            if level == 1 or level == 2:
                while len(fail_list) > (27*(1-prob_survive)):

                    # Pick a random box where no star formed
                    random_box = np.random.randint(0, len(fail_list))

                    # Get the positions of the random box 
                    dx = fail_list[random_box][0]
                    dy = fail_list[random_box][1]
                    dz = fail_list[random_box][2]

                    # Try again to create a star there
                    success, new_box_positions, new_box_vels = get_pos_and_vel(box_positions, box_vels, box, new_box_positions, new_box_vels, prob_survive, side_length, level, dx, dy, dz)

                    # If a star is made there remove this position's coordinates from
                    # the list of boxes where stars failed to form
                    if success:
                        fail_list.pop(random_box)
   
        side_length = side_length / 3.
        
        box_positions = new_box_positions
        box_vels = new_box_vels


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #cut sphere
    for star in range(len(box_vels[0])-1, -1, -1):

        if np.sqrt(box_positions[0][star]**2 + box_positions[1][star]**2 + box_positions[2][star]**2) > 0.8:

            box_positions[0].pop(star)
            box_positions[1].pop(star)
            box_positions[2].pop(star)
            box_vels[0].pop(star)
            box_vels[1].pop(star)
            box_vels[2].pop(star)


    #Cull excess stars
    #calc how many stars to delete

    if len(box_vels[0]) < n_stars:
        print('too little')
        desired_number = len(box_vels[0])
    else:
        desired_number = n_stars
    to_delete = len(box_vels[0]) - desired_number

    #Choose random stars to delete
    remove_list = range(len(box_vels[0]))
    remove_list = np.random.choice(remove_list, size=to_delete, replace=False)
    remove_list = sorted(remove_list, reverse=True)

    for star in remove_list:
        box_positions[0].pop(star)
        box_positions[1].pop(star)
        box_positions[2].pop(star)
        box_vels[0].pop(star)
        box_vels[1].pop(star)
        box_vels[2].pop(star)


    #Add a bit more noise to star positions and velocitiies
    box_positions[0] = [np.random.normal(loc=pos, scale=0.05) for pos in box_positions[0]]
    box_positions[1] = [np.random.normal(loc=pos, scale=0.05) for pos in box_positions[1]]
    box_positions[2] = [np.random.normal(loc=pos, scale=0.05) for pos in box_positions[2]]
    box_vels[0] = [np.random.normal(loc=vel, scale=0.05) for vel in box_vels[0]]
    box_vels[1] = [np.random.normal(loc=vel, scale=0.05) for vel in box_vels[1]]
    box_vels[2] = [np.random.normal(loc=vel, scale=0.05) for vel in box_vels[2]]
    
    # If velocity structure is false then shuffle the velocities
    # in a random order to remove any structure
    if not vel_struct:
        rand_order = [i for i in range(n_stars)]
        np.random.shuffle(rand_order)
        box_vels[0] = [box_vels[0][i] for i in rand_order]
        box_vels[1] = [box_vels[1][i] for i in rand_order]
        box_vels[2] = [box_vels[2][i] for i in rand_order]

    # Convert to numpy arrays
    r = np.array(box_positions)
    v = np.array(box_vels)

    return r, v

def get_masses(n_stars):

    # Create an array to hold them
    m = np.zeros((n_stars))
    
    # Get masses for each star, drawn from the Maschburger IMF.
    m = maschberger_IMF.maschberger_IMF(n_stars, 0.1, 30)

    return m


# SET UP VARIABLES AND ARRAYS NEEDED TO MAKE THE FRACTAL

# Read the input values from the command line
try:
    n_stars = int(sys.argv[1]) 
    prob_survive = float(sys.argv[2]) 
    desired_half_mass_radius = float(sys.argv[3]) 
    
except:
    
    print('Error: Need to specify the number of stars, required probability of survival and desired half mass radius via the command line.')
    print('For example a valid command would be: python2 make_box_fractal.py 1000 0.5 5')
    
    exit()

# Check the probability of survival is between 0 and 1. If not throw an error
if (prob_survive <= 0) or (prob_survive > 1):
    
    print('Error: the probability a star survives to the next fractal generation must be between 0 and 1')
    exit() 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# BUILD THE FRACTAL

# Make the fractal to get positions and masses for the stars 
r, v = make_fractal(n_stars, prob_survive)

# Get masses for the stars
m = get_masses(n_stars)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SCALE THE POSITIONS TO THE DESIRED RADIUS

# First need to move to the centre of mass
# Move the centre of mass to zero
r = CoM(n_stars, r, m)

# Find the half mass radius of the cluster
half_mass = sum(m) / 2.
half_mass_radius = 0.
mass_inside = 0.

#Keep increasing half_mass_radius until the mass inside exceeds the half mass
while (mass_inside < half_mass):

    half_mass_radius += 0.05
    mass_inside = 0.

    for star in range(n_stars):

        dist_from_center = np.sqrt(r[0][star]**2 + r[1][star]**2 + r[2][star]**2)
        if dist_from_center < half_mass_radius:

            mass_inside += m[star]


# Scale all the positions so the half mass radius is the desired one
r = r * (desired_half_mass_radius / half_mass_radius)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE THE FRACTAL TO FILE

# Transpose the arrary so the x y z data is saved as 3 columns
# Rather than 3 long rows.
r = r.T
v = v.T

# Save to text files
np.savetxt('r_stars.txt', r)
np.savetxt('v_stars.txt', v)
np.savetxt('m_stars.txt', m)

