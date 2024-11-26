# RuO2_OER_MLP
The repository for MLP train and OER investigation.

**2024.10.25**
- Modify the candidates.
- Modify the generations
- How to link the calculator to generations.
- How to set learning workflow.

**2024.10.26**
- Modify the candidates.
    - When pop some atom in dissolution, the fixed atom will distrub.
    - get_surface() will get stuck in loop due to find the metal atoms.
- Growth method can not success in vacuum. It need a atom, it is not correct.
- Grow is tend to attach to the new grow atom, due to the get_surface is depend on the z-zxis.
