# attitude

C++ library containing lightweight implementations of various attitude description sets.

It is assumed that whomever may use this library does have a knowledge of attitude representation, the mathematics involved and the relevant trade-offs. A cheatsheet detailing this to some degree is to come (WIP). This is intended to provide a method of working with attitude parameterisations without cluttering higher level programs with the specifics of the kinematics.

Currently implemented:  
 - Euler Angles (any order)  
 - Quaternions (or Euler Parameters)  
 - Classical Rodriguez Parameters
 - Modified Rodriguez Parameters

_(Of course, rotations can also be expressed as directional cosine matrices and principal rotation vectors. Rudimentary linear algebra library included.)_

Each representation has the following encoded:
 - Addition of rotations can be completed simply using `operator+,+=`. Relative rotaions can be computed with `operator-,-=`. Linear algebra operations (direct where possible - otherwise via the rotational matrices) are encoded.
 - Propagation in time: each implementation can be propagated across a timestep if the angular veloctiy is known.  
 - Easy conversions from one to another.

This library has been written with the following considerations:
 - Header-only: less administrative overhead to include in larger projects. 
 - Templated: Precision of parameter sets is free to be specified based on performance limitation/desired accuracy.  
 - Lightweight: avoided use of Boost/STL where possible, memory allocation at compilation.
