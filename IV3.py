import os

# values of the applied oltages
drain_voltages=[0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
gate_voltage=[0.0, 0.5, 1.0, 1.5, 2.0]

# for every applied voltages, create the script for Archimedes
# and run import
for gate in gate_voltage:
    for drain in drain_voltages:
        f = open("script.input","w")
        script = "\
        TRANSPORT MC ELECTRONS\n\
        MOSFET\n\
        FINALTIME 15.0e-12\n\
        TIMESTEP 0.0015e-12\n\
        XLENGTH 0.2e-6\n\
        YLENGTH 3.0e-9\n\
        XSPATIALSTEP 120\n\
        YSPATIALSTEP 30\n\
        ACOUSTICSCATTERING ON\n\
        OPTICALSCATTERING ON\n\
        MATERIAL X 0.0 0.2e-6    Y 0.0 3.0e-9  GERMANIUM\n\
        CONDUCTIONBAND PARABOLIC\n\
        DONORDENSITY    0.       0.         0.2e-6    3.0e-9   1.e21\n\
        DONORDENSITY    0.       2.25e-9    0.034e-6    3.0e-9   5.e21\n\
        DONORDENSITY    0.166e-6   2.25e-9    0.2e-6    3.0e-9   5.e21\n\
        ACCEPTORDENSITY 0.       0.         0.2e-6    3.0e-9   1.e21\n\
        CONTACT DOWN  0.0    0.2e-6 INSULATOR 0.0\n\
        CONTACT LEFT  0.0    3.0e-9 INSULATOR 0.0\n\
        CONTACT RIGHT 0.0    3.0e-9 INSULATOR 0.0\n\
        CONTACT UP    0.034e-6 0.068e-6 INSULATOR 0.0\n\
        CONTACT UP    0.132e-6 0.166e-6 INSULATOR 0.0\n\
        CONTACT UP    0.0    0.034e-6 OHMIC     0.0 5.e21\n\
        CONTACT UP    0.068e-6 0.132e-6 SCHOTTKY  "+str(gate)+"\n\
        CONTACT UP    0.166e-6 0.2e-6 OHMIC       "+str(drain)+"  5.e21\n\
        QEP_MODEL DENSITY_GRADIENT\n\
        LATTICETEMPERATURE 300.\n\
        STATISTICALWEIGHT 200\n\
        MEDIA 3\n"
        f.write(script)
        f.close()
        command="./a.out script.input > output_"+str(gate)+"_"+str(drain)+".txt"
        os.system(command)

# remove the script
os.system("rm -rf script.input")