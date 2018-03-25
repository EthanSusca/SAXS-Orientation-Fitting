# SAXS-Orientation-Fitting
Finds Bragg peaks from 2D SAXS .tiff file and determines best-fit Euler Angles (Bunge, passive rotations).

General Structure

WriteEBSDtxt: Iterates through all 2D SAXS image files to generate text file that can be interpreted by EBSD analysis software

  Euler Image:  With the parameter set in WriteEBSDtxt, this script first iteratively identifies Bragg peak locations and then fits them
    Fit_Diffraction_Spot  (see below)
    Fit_Orientation_Bunge (see below)

Fit_Diffraction_Spot:       Bragg Peak Location Fitting   
For a given sample, high intensity pixels (approximate centers of Bragg reflections) were identified, fit with a polar two-dimensional Gaussian (a function of the radial and azimuthal directions) to determine the center of the reflection. The centers of reflections were recorded for all reflections with intensity larger than a pre-determined threshold value (set in Euler Image). Two separate threshold intensities were employed: one for low q reflections (which include the {211} and {220} Bragg peaks) because the background at low q is on the order of 100–1000, and another threshold for all high q Bragg reflections (which had a maximum background with intensity of 10–20). After each reflection location was recorded, all pixels +/- two standard deviations from center were assigned a value of zero. The Bragg reflection fitting process was iteratively employed until no pixels were identified with intensity larger than the lower threshold value.

Fit_Orientation_Bunge:  	Orientation Fitting, Euler Angle Identification
A small angle approximation, where the Ewald sphere is assumed to be a plane, was used to calculate expected Bragg peak reflections for a given set of Euler angles using the Bunge convention with passive rotations. An orientation space spanning 90˚ x 90˚ x 90˚ was explored, corresponding to 6 symmetrically equivalent orientation spaces. A mosaicity parameter was set to 8˚ so that all Bragg reflections within 8˚ of the Ewald plane were established as calculated reflections for a specific orientation. Experimentally identified Bragg reflection locations within 5% of the calculated Bragg peak locations (with respect to q in radial and azimuthal directions) were recorded as a successful fit of an experimental reflection to a calculated reflection. The set of Euler angles that yielded the greatest number of agreements between calculated and experimental reflections were used to identify the best-fitting orientation and the Euler angles associated with this orientation were recorded. If two orientations had the same number of experimental/calculated agreements, the orientation with the smallest residual was recorded. The set of Euler angles, number of successfully fit reflections, and the fraction of experimental reflections accounted for by the fit were recorded for all 2D SAXS pattern acquired. The commented code used to analyze the data is included in the supporting text at the end of this document.
