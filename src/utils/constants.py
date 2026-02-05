# constants.py

MSUN_TO_G = 1.989e33       # Solar Mass to grams
AU_TO_CM = 1.496e13        # AU to centimeters
PC_TO_CM = 3.0856e18       # Pc to centimeters
W_TO_ERG_S = 1e7           # Watts to ergs/s

# (Msun / AU^2) -> (g / cm^2)
MSUN_AU2_TO_G_CM2 = MSUN_TO_G / (AU_TO_CM ** 2)
# (Msun / PC^2) -> (g / cm^2)
MSUN_PC2_TO_G_CM2 = MSUN_TO_G / (PC_TO_CM ** 2)
# (Msun / AU^3) -> (g / cm^3)
MSUN_AU3_TO_G_CM3 = MSUN_TO_G / (AU_TO_CM ** 3)
# (Msun / PC^3) -> (g / cm^3)
MSUN_PC3_TO_G_CM3 = MSUN_TO_G / (PC_TO_CM ** 3)
# AU -> PC
AU_TO_PC = AU_TO_CM / PC_TO_CM

# Jansky to CGS specific flux density (erg/s/cm^2/Hz)
JY_TO_CGS = 1e-23