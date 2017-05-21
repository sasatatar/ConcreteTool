# ConcreteTool

### MATLAB program for reinforced concrete section design

This program was created as a result of a [diploma thesis](https://drive.google.com/open?id=0Byvv9xAhq3xaLVJ2Wi1LNVpxSTA) (in Serbian).

ConcreteTool can be used to design reinforced concrete beam sections (rectangular and T sections) subjected to bending, shear and torsion loads, according to Eurocode 2. The program has easy to use graphical user interface.

The program requires at least MATLAB R2014b, as well as polyxpoly from Mapping Toolbox (http://www.mathworks.com/help/map/ref/polyxpoly.html).
If you don't have Mapping Toolbox installed, you can use this:
https://gist.github.com/ashaindlin/75e71b2825bb112fdbf9

To start the program GUI, save all files in one folder (keeping the subfolder `section_templates`) and set it as current directory in MATLAB. 

You can then start the program by typing the command `ConcreteTool` inside Command window in Matlab.

### Program sceenshots

#### Define cross-section

![image](https://cloud.githubusercontent.com/assets/5138412/26280373/1dc166ea-3dd1-11e7-868f-ae9a26392d3d.png)

#### Set material properties and loads

![image](https://cloud.githubusercontent.com/assets/5138412/26280362/bf23b34a-3dd0-11e7-9653-822ec4d4bec7.png)

#### Over reinforced cross-seciton

![image](https://cloud.githubusercontent.com/assets/5138412/26280377/427e1488-3dd1-11e7-8b35-7ef403921dc6.png)

#### Add/remove rebars

![image](https://cloud.githubusercontent.com/assets/5138412/26280375/34f2018a-3dd1-11e7-84c3-a166541cb3c0.png)

#### M-phi (curvature) diagram

![image](https://cloud.githubusercontent.com/assets/5138412/26280379/4986336e-3dd1-11e7-90c0-71bd464e597a.png)

#### Vrd section design

![image](https://cloud.githubusercontent.com/assets/5138412/26280384/53ec3a7e-3dd1-11e7-9dfe-af25233424df.png)

#### Vrd + Trd section design

![image](https://cloud.githubusercontent.com/assets/5138412/26280389/8141bfbc-3dd1-11e7-8c68-72382c8ea6dc.png)

#### Vrd - Trd interaction diagram

![image](https://cloud.githubusercontent.com/assets/5138412/26283084/1f2038b8-3e21-11e7-8d0e-e0dff5dda3b9.png)
