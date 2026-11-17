"""Training on real antibody-antigen structures.

dataset.py : structures -> cropped/centered/padded tensors (CDR design block +
             antigen context), with a generate_mask for inpainting.
train.py   : the training loop (flow matching + CTMC) for the IPA CoDesignModel,
             with validation and checkpointing.
"""
