import BeamModule

beamFromAperture = BeamModule.Beam.initialiseBeamFromApMeas("20141216/20141216_Beam_profile_x_sma_ap.dat","20141216/20141216_Beam_profile_y_sma_ap.dat","threshold", 10, 2, "beamtest.pgm")
beamFromPNGImage = BeamModule.Beam.initialiseBeamFromPNG("beam_z20_g1_e70_s100_60_t100.png",0.3,0.6,0.1,"anotherbeamtest.pgm")
beamFromPGM = BeamModule.Beam.initialiseBeamFromPGM("jonny.pgm")
