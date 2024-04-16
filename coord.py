#%%
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import pandas as pd

# 主程序

def convCoord(coord, frame_in, frame_out=None):
    if frame_in == "J2000":
        frame_in = "fk5"
    elif frame_in == "gal":
        frame_in = "galactic"
    
    c = SkyCoord(coord, frame=frame_in, unit=(u.hourangle, u.deg) if ":" in coord else (u.deg, u.deg))
    
    if frame_out is None:
        c_out_gal_deg = f"galactic in deg: {c.transform_to('galactic').l.deg:.2f}, {c.transform_to('galactic').b.deg:.2f}"
        c_out_fk5_hms = f"fk5 in hms: {c.to_string('hmsdms').split()[0]}, {c.to_string('hmsdms').split()[1]}"
        c_out_fk5_deg = f"fk5 in deg: {c.transform_to('fk5').ra.deg:.2f}, {c.transform_to('fk5').dec.deg:.2f}"
        return '\n'.join([c_out_gal_deg, c_out_fk5_hms, c_out_fk5_deg])
    elif frame_out == "fk5_hms":
        return f"fk5 in hms: {c.to_string('hmsdms').split()[0]}, {c.to_string('hmsdms').split()[1]}"
    elif frame_out == "fk5_deg":
        return f"fk5 in deg: {c.transform_to('fk5').ra.deg:.2f}, {c.transform_to('fk5').dec.deg:.2f}"
    elif frame_out == "galactic":
        c_out_fk5_deg = f"fk5 in deg: {c.transform_to('fk5').ra.deg:.2f}, {c.transform_to('fk5').dec.deg:.2f}"
        return c_out_fk5_deg
    else:
        return "Something wrong about coord frame"

def convCoord_from_list(coordFile, frame_out="fk5"):
    df = pd.read_csv(coordFile, header=None, names=["c1", "c2", "frame"], sep=" ")
    for index, row in df.iterrows():
        print(f"We are working on {index+1}")
        c1 = row["c1"]
        c2 = row["c2"]
        frame_in = row["frame"]
        coord = f"{c1} {c2}"
        out = convCoord(coord, frame_in, frame_out)
        print(out)

coordFile = "./coord_list.txt"
convCoord_from_list(coordFile)
