=================
Git quick summary
=================


To setup git repo on a new machine, clone the bitbucket repo :: 
    
    $ git clone git@bitbucket.org:py_gem/pyclusterevol.git

Now tell git who you are, and that you have an access to the repository. ::

    $ git init... I forgot! :)

Then to/from which branch you will push/pull. ::

    $ git --set-upstream remote_name branch_name


To discard modification to a file/files ::

    $ git checkout -- file.name(s)

For example, if you made a change to files on a server and pushed the modification.
But forgotten doing so, you again start to modify a file on your laptop and realize that you are doing something redundant and you could simply pull from the repository. Then you can use "chekcout --" to revert the file so that you can pull from the repo. 



