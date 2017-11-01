# Filament
Program for pair-wise superposition of neighbouring domains in multidomain proteins. For using it you need to have ATSAS package (https://www.embl-hamburg.de/biosaxs/download.html) installed on your computer.
The program performs structural alignment of domains and provides one with tranformational matrix (as well as Euler angles) for each pair of domains. User have two possibilities how to perform the domains subdivision: using automatic regime based on NMA (PARCOOR program), or manual indication of borders between different modules. The superposition is performed by SUPCOMB program using NSD metrics. Additionally, Filament calculates Twist and Tilt angles according to [Peer Bork et al, 1996].
At the moment it works only under Unix OS.
