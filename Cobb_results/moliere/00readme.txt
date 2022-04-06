
18/V/17

The file names for the Moliere calculation are of the form:

  Material-thickness-momentum-mm-rf.dat

where:

   thickness = thickness in cm (so 6p5 means 6.5 cm)

   momentum = p_\mu in MeV/c

   mm = Mol = space angle
      = Molprj = projected angle

   rf = raw = at the points given by Moliere
      = fine = using my cubic spline interpolation

The files are plain ascii (text) files where:

   column 1 is angle (radians)

   column 2 is P, probability per radian at this angle.

The 'fine' files have 800 points (bins):

    0.0 to 0.2 radians for the space angle distributions
   -0.2 to 0.2 radians for the projected angle distributions

To normalise to the data, P --> P * total number of events * bin size

where the total number of events is the number of _upstream_ muons (and the downstream data corrected for efficiencies etc.).

P is calculated at the bin centres for the 'fine' (interpolated) files.

The 'raw' distributions are simply calculated at the angles given by Moliere.

For the record:

 I used the 3D tables of Bethe and Moliere's tables for the projected angles.

 I used the Fano correction for electrons as recommended by Gotshchalk with -u = -5.0.

The file example.pdf shows what the 'fine' distributions look like for 35 cm of Xenon with the raw points superimposed.



