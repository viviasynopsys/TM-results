# lumerical/

The three Lumerical project files saved by the mixed-grating run (`rcwa_mixed_A.fsp`, `rcwa_mixed_B.fsp`, `rcwa_mixed_C.fsp`, about 27 MB together) are **not included in this commit**. They are saved copies of Ansys's "Volume Holographic Grating" example model, which is copyright Ansys, Inc. and subject to the Ansys Terms of Service (see `remote_pack_conical/lumerical/ansys_example/license.txt`).

They are kept in the local `F-PSM-mixed/lumerical/` folder; the driver writes `lumerical/rcwa_mixed_<set>.fsp` for the first config it solves in each set (`python lumerical_rcwa_mixed.py --set A|B|C`). Everything else the run produced is here.
