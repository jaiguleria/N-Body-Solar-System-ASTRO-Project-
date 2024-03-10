#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from astropy.time import Time
from astropy.coordinates import get_body_barycentric, get_body, get_moon
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric_posvel
solar_system_ephemeris.set('jpl')

#Setting an initial time, an empty array to store initial positions and the chosen planets to run the loop for
t = Time("2000-09-22 12:00")
initialpos = []
planets = ["sun", "earth-moon-barycenter", "mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune"]

#Working on my for loop that will calculate the positions and velocities using AstroPy
for planet in planets:
    pos = get_body_barycentric_posvel(planet, t)
    ggupos = [pos[0].x.to("AU").value,pos[0].y.to("AU").value,pos[0].z.to("AU").value,pos[1].x.to("AU/d").value,pos[1].y.to("AU/d").value,pos[1].z.to("AU/d").value]
    initialpos.append(ggupos)
    print(f"{planet}: x={ggupos[0]}, y={ggupos[1]}, z={ggupos[2]}, dx={ggupos[3]}, dy={ggupos[4]}, dz={ggupos[5]}")
    

#Get the barycentric coordinates of the solar system
barycenter_coords = get_body_barycentric('sun', t)

print("Barycenter coordinates (x, y, z) in AU:", barycenter_coords)


# In[3]:


#Defining my n-body function
def N_body_function(t, x, masses, ai):

    dXN = np.zeros(6*len(masses))  #initialize acceleration array
    c = 3*(10**8)
    
    for i in range(len(masses)): #first body
        
        dXN[6*i:(6*i)+3] = x[(6*i)+3:(6*i)+6] #set velocities as first derivatives of positions
        
        for j in range(len(masses)): #second body
            
            if j != i: #if second body isn't the same as the first one
              
                rn = np.sqrt((x[6*j]-x[6*i])**2+(x[(6*j)+1]-x[(6*i)+1])**2+(x[(6*j)+2]-x[(6*i)+2])**2) #calculate the separation between bodies
                
                if i != 0 and j == 0: #if planet isn't sun
                
                    force = - masses[j] * ((x[6*j:6*j+3] - x[6*i:6*i+3])/rn**3)*(1 - (9*GM/(c**2*ai[i]) + 6*GM/c**2*rn)) #calculate gravitational force exerted by body j on body i, with the GR correction applied
                
                else:
                    
                    force = - masses[j] * ((x[6*j:6*j+3] - x[6*i:6*i+3])/rn**3) #for the sun, no GR correction
                
                dXN[6*i+3:6*i+6] += force #update acceleration components of body i, += updates instead of overwriting
    
    return dXN


#Planet mass ratios
sem = 328900.56
mer = 6023600
ven = 408523.71
mar = 3098708
jup = 1047.3486
sat = 3497.898
ura = 22902.98
nep = 19412.24
GM = -(0.01720209895)**2


#Setting up my initial position & velocity array, alongisde my GM/mass values array
initialpos = np.array(initialpos).flatten()
print(initialpos)
N_Mass = np.array([1, 1/sem, 1/mer, 1/ven, 1/mar, 1/jup, 1/sat, 1/ura, 1/nep])*GM
ai = [0, 0.38700, 0.72300, 1, 1.52400, 5.20440, 9.58260, 19.21840, 30.11000] #semi-major axes


#Call solve ivp and set time span
times = [0,7300.0] 
N_Body = solve_ivp(N_body_function, times, initialpos, args=(N_Mass, ai),method='RK45',rtol=1e-12) 
from mpl_toolkits.mplot3d import Axes3D

#Plotting our results & establishing colour & label arrays
fig = plt.figure(figsize=(14,12))
ax = fig.add_subplot(111, projection='3d')
colours = ['y', 'b', 'gray', 'orange', 'r', 'brown', 'pink', 'green', 'black']
labels = ["Sun", "Earth", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]

#Automating the plotting process
for i in np.arange(len(planets)):
    ax.plot(N_Body.y[6*i], N_Body.y[6*i+1], N_Body.y[6*i+2], color = colours[i], label = labels[i]) 

ax.set_xlabel('AU (m)')
ax.set_ylabel('AU (m)')
ax.set_zlabel('AU (m)')
ax.legend()
plt.title("N-Body System")


# ## Session 7 - Finish GR, JPL Comparisons and begin final project

# To start off with, I will refamiliarise myself with the code for the GR correction and attempt to finish this and half of the JPL comparison by lunchtime. 
# 
# - Looking ahead to JPL, I'm setting a time span of 20 years of data stretching from 2000-09-22 to 2020-09-22, this seems like a reasonable amount of time in which to record results as 
# - Figured out that for my GR correction, my semi-major axes array was not even big enough to match up to all the planetary masses and since the Sun won't have a GR correction for itself, its semi-major axis will be 0. In my code: `if i != 0: force = - masses[j] * ((x[6*j:6*j+3] - x[6*i:6*i+3])/rn**3)*(1 - (9*GM/(c**2*ai[i]) + 6*GM/c**2*rn))` I need to add an `else` to account for the sun having no GR correction ie. an index of 0.
# - I also ran into trouble with the JPL horizons ephemeris generator as I tried to generate it with a geocentric system, whereas I should be doing it with the solar system barycentre. This isn't even necessary for the comparison, as I can simply use the `astropy planetary positions` file's code instead.
# - Thought I fixed the GR correction but I ran into another problem with the straight lines! I plan to investigate the gravitational aspect of it as the minus sign relating to gravity is what took up my entire session last week and this seems to be a similar error. Turns out ` else: force = - masses[j] * ((x[6*j:6*j+3] -x[6*i:6*i+3])/rn**3) dXN[6*i+3:6*i+6] += force #update acceleration components of body i, += updates instead of overwriting` was just updating the force array only after the code written specifically for the Sun, this was fixed simply by removing the indent.
# - I'm looking ahead to the JPL comparisons for now, seeing the difference between the `astropy planetary positions` code taking the positions straight from the solar system ephemeris, against my own planets code. Both start at the same initial conditions, it's just that the JPL code will generate a list of times/positions that will be slightly different from my N-body, I need to produce a list of these times/positions and plot these against my code. `Solve ivp` gives times from just zero instead of from the date we set like in `astropy` code. So what I wanna do is add the time from the date and time set in my `t = Time("2000-09-22 12:00")` but converted to `Jd`, to my solver which starts from 0 in Jd.
# - `i` is planet, `0, 1, 2` take x, y, z and all the times stretching to the right of the array

# In[32]:


#![fail.PNG](fail.PNG) INSERT IMAGE HERE IN REPORT AS EXAMPLE OF A HURDLE


# In[25]:


#setting up our time in julian days
tJD = Time("2000-09-22 12:00")
timejd = tJD.jd + np.array(N_Body.t)

#making them into JD time objects instead of string/integers
TIME = Time(timejd, format='jd', scale='utc') 

posvel_bodies = []

for planet in planets: #pos defined within loop so no global variable issues
    
    pos = get_body_barycentric_posvel(planet, TIME)
    ggupos = [pos[0].x.to("AU").value,pos[0].y.to("AU").value,pos[0].z.to("AU").value,pos[1].x.to("AU/d").value,pos[1].y.to("AU/d").value,pos[1].z.to("AU/d").value]
    posvel_bodies.append(ggupos)


# In[14]:


#plotting our JPL ephemeris data in the same way as our solver plot
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111, projection ='3d')

for i in np.arange(len(planets)):
    ax.plot(posvel_bodies[i][0], posvel_bodies[i][1], posvel_bodies[i][2], color = colours[i], label = labels[i]) #no .y as that's only for solvers and we have a 3d array

ax.set_xlabel('AU (m)')
ax.set_ylabel('AU (m)')
ax.set_zlabel('AU (m)')
ax.legend()
plt.title("JPL Ephemeris")


# - I am on track, I now want to plot the differences between the datasets on some graphs and then speak to a demonstrator about my project in which I want to plot Voyager's orbit and quantify a result in some way
# - Need to plot the differences as an `r` calculation between the two bodies at every point, in a `for` loop and adding the x, y, z coordinates together to get a single value that I can plot against the time
# - Trying to figure out how to let the code differentiate between my solver's data and the ephemeris JPL data with the use of indices yet again. I'm now trying to deduce the shape of the `posvel_bodies` array so that I can plot my data,
# - Need to plot it against time retrieved from the N-Body solver itself, also need to set up an array to store all these datapoints in and make sure to not overwrite what I'm calculating each time.
# - After some slight trouble with subplotting, I have figured out how to plot the differences.

# In[120]:


#plotting the data
diff_figure = plt.figure(figsize = (25, 35))

#set up an empty list to store differences
diff_array = []  

#setting up my for loop to calculate the differences
for i in range(len(planets)):
        
    difference = np.sqrt((N_Body.y[6*i] - posvel_bodies[i][0])**2 + (N_Body.y[(6*i)+1] - posvel_bodies[i][1])**2 + (N_Body.y[(6*i)+2] - posvel_bodies[i][2])**2)
        
    diff_array.append(difference)
    
    #add a subplot for the current planet
    ax = diff_figure.add_subplot(len(planets), 3, i+1)
    
    #plot the difference on the current subplot
    ax.plot(N_Body.t, difference, color = colours[i], label=labels[i])
    ax.set_ylabel('Position Difference (AU)')
    plt.xlabel('Time from [tJD] (days)')
    ax.set_title(f'Difference in position for {labels[i]}')
    plt.grid()

#set common x-axis label and title
plt.suptitle('Difference in position (AU) between our Solver and JPL Ephemeris Data'
             '\n'
             '\n', fontsize=25)
diff_figure.subplots_adjust(top=0.7)
plt.tight_layout()
plt.show()

