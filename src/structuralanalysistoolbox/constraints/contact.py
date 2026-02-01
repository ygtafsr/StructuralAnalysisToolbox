
"""
3.2. Steps in a Contact Analysis
Following are basic steps for performing a typical surface-to-surface contact analysis:
1. Create the model geometry and mesh (p. 16)
2. Identify the contact pairs (p. 17)
3. Designate contact and target surfaces (p. 18)
4. Define the target surface (p. 20)
5. Define the contact surface (p. 27)
6. Set the element KEYOPTS and real constants (p. 35)
7. Define/control the motion of the target surface (rigid-to-flexible only) (p. 110)
8. Apply necessary boundary conditions (p. 111)
9. Apply fluid pressure-penetration loads (p. 112)
10. Define solution options and load steps (p. 118)
11. Solve the contact problem (p. 127)
12. Review the results (p. 128)

Target and contact elements that make up a contact
pair are associated with each other via a shared real constant set.

Different contact pairs must be defined by a different real constant set, even if
the element real constant values do not change.
"""

class Contact:
    pass