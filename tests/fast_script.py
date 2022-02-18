from rimopy import spice_iface as spi
KERNELS_PATH = "./kernels"
VALL_LAT = 41.6636
VALL_LON = -4.70583
VALL_ALT = 705
JAN_FULL_MOON_00 = "2022-01-17 00:00:00"
JAN_FULL_MOON_17 = "2022-01-17 17:00:00"
FEB_NEW_MOON_00 = "2022-02-02 00:00:00"
md = spi.get_moon_data(VALL_LAT, VALL_LON, VALL_ALT, JAN_FULL_MOON_00, KERNELS_PATH)
print(md.absolute_mpa_degrees)
