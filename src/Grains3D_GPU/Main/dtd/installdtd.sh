#!/bin/bash
# Generates the built DTD (Grains3D.dtd) from the source template by substituting
# absolute paths for the Convex and Contact entity placeholders.
# Requires GRAINS_HOME to be set (source Env/grainsGPU.env.sh first).

cp Generic_Grains3D.dtd Grains3D.dtd

sed -i "s|___WHERESISCONVEX__|${GRAINS_HOME}/Main/dtd/Convex.dtd|g"  Grains3D.dtd
sed -i "s|___WHERESISCONTACT__|${GRAINS_HOME}/Main/dtd/Contact.dtd|g" Grains3D.dtd
