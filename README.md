# elisetools

A repo to store helpful tools I've made for fMRI analyses with the CanLabCore tools and SPM (see links below). 
I've been making what is useful for me, so feel free to edit for your own needs. 

## Setup

Add to your `startup.m`:

```matlab
addpath(genpath('/path/to/elisetools'))
```

## Tools - just one for now but maybe more added soon who knows

### neuroimaging/check_mask_coverage 
- When applying masks for pattern expression, it is easy to inadvertently use a mask that was defined in a different MNI space, at a different resolution, or with a different field of view than your data. This function checks if a binary mask adequately covers a beta/contrast image. 
- It can tell you if its too small, or too big, with options for plotting the overlap (or lackthereof).
- You can also skip plotting and just compare the headers of two .nii files (e.g. comparing data with a collaborator).
- It was coded using Claude and tested on my data, but I tried to make helpful edits to the text and code to make it more readable for someone learning neuroimaging. 

To do: add image of the output 

## Requirements

- MATLAB
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/)
- [CANlab Core Tools](https://github.com/canlab/CanlabCore)

## License

MIT License — see [LICENSE](LICENSE) for details.
