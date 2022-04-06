# Workflow analysis

Still in progress ...

## Descriptor analysis

To run the descriptor analysis just update the `descriptor.toml` file with the desired parameters.

Then run:

`python frame_by_frame_descriptor.py -c descriptor.toml` 

to get the descriptor output in a `.npy` compressed file.


## Low dimensional embedding

To run the code for reducing the dimensionality of a dataset update the `embed.toml` file with the desired parameters.

Then run:

`python lowDim_embedding.py -c embed.toml` 

to get the data embedded in a low dimensional manifold output in a `.npy` compressed file.