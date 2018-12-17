import os

# values of the applied oltages
applied_voltages=[0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]

# for every applied voltages, create the script for Archimedes
# and run import
for voltage in applied_voltages:
    f = open("script.input","w")
    script = "\
    TRANSPORT MC ELECTRONS\n\
    FINALTIME 15.0e-12\n\
    TIMESTEP 0.0015e-12\n\
    XLENGTH 0.2e-6\n\
    YLENGTH 3.0e-9\n\
    XSPATIALSTEP 120n\
    YSPATIALSTEP 120\n\
    ACOUSTICSCATTERING ON\n\
    OPTICALSCATTERING ON\n\
    MATERIAL X 0.0 0.2e-6    Y 0.0 3.0e-9  GERMANIUM\n\
    CONDUCTIONBAND PARABOLIC\n\
    DONORDENSITY    0.       0.         0.2e-6    3.0e-9   1.e22\n\
    DONORDENSITY    0.       2.25e-9    0.034e-6    3.0e-9   5.e22\n\
    DONORDENSITY    0.166e-6   2.25e-9    0.2e-6    3.0e-9   5.e22\n\
    ACCEPTORDENSITY 0.       0.         0.2e-6    3.0e-9   1.e19\n\
    CONTACT DOWN  0.0    0.2e-6 INSULATOR 0.0\n\
    CONTACT LEFT  0.0    3.0e-9 INSULATOR 0.0\n\
    CONTACT RIGHT 0.0    3.0e-9 INSULATOR 0.0\n\
    CONTACT UP    0.034e-6 0.068e-6 INSULATOR 0.0\n\
    CONTACT UP    0.132e-6 0.166e-6 INSULATOR 0.0\n\
    CONTACT UP    0.0    0.034e-6 OHMIC     0.0 5.e22\n\
    CONTACT UP    0.068e-6 0.132e-6 SCHOTTKY  -1.0\n\
    CONTACT UP    0.166e-6 0.2e-6 OHMIC       "+str(voltage)+"  5.e22\n\
    OXYDE   UP    0.068e-6 0.132e-6 0.005e-6 0.3\n\
    OXYDE   DOWN  0.0      0.2e-6  0.010e-6 0.3\n\
    QEP_MODEL DENSITY_GRADIENT\n\
    MAXIMINI\n\
    LATTICETEMPERATURE 300.\n\
    STATISTICALWEIGHT 200\n\
    MEDIA 3\n"
    f.write(script)
    f.close()
    command="./a.out script.input > output_"+str(voltage)+".txt"
    os.system(command)

# remove the script
os.system("rm -rf script.input")