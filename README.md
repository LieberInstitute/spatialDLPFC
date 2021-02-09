# spython

Testing ground (sandbox) for `libd_jhu_spatial_python`.

# Notes on organization & communication

* None of the code in this repo will be made public (as is). This repo is solely for testing external (python) tools.
* Use full paths in your scripts instead of relative paths. Doing so will enable us to move our test code from this repo to the repos used for the analyses of our projects.
* Try out methods with their example data first: this involves installing the software and adjusting any configuration required for JHPCE (like making a JHPCE module). If it works, then try it with a subset of samples from one of our spatial projects (some methods might work better with the data from a brain region than other regions given the different anatomical properties). If we like it, then we'll move the code to the real analysis repo and apply it there with the rest of the samples for that project.
* Keep track of the order of the commands you run. Some `conda` commands at JHPCE reset environment variables, which can play a crucial role with python dependencies, and thus for understanding bugs.
* Make git commits often.
* Use github issues to keep track of our progress and any related external interactions (issues on other repos, emails, etc).
* When interacting with others, issues on github repos are preferred over emails. That way others can "see" the interactions and chime in. Similarly, we prefer that you use the `libd_jhu_spatial_python` Slack channel over direct messages.
* Prior to creating an issue on an external repo, it's best if you draft the issue and meet to discuss the draft, such that you can check it for clarity and completeness. Sometimes when you write an issue, you find the solution to the problem yourself.

# JHPCE

Location: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spython`

* Use the `lieber_jaffe` user group.
* Remember to edit your `~/.bashrc` file as shown at https://lcolladotor.github.io/bioc_team_ds/config-files.html#bashrc to include the command `umask u=rwx,g=rwx,o=`. This will make our lives easier in terms of working with the same JHPCE clone.
