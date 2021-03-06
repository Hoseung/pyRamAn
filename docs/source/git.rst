=================
Git quick summary
=================

Set up a new repo
-----------------
*
If you are setting up Git on a new machine, 
first you need to register the machine to your bitbucket profile.
https://confluence.atlassian.com/bitbucket/set-up-ssh-for-git-728138079.html,

0. set up a ssh key set on the local machine and copy the key.::
    $ cd ~/.ssh
    $ ssh-keygen
    $ cat ~/.ssh/id_rsa | xclip  

1. add the key to (bitbucket avata) - Bitbucket settings - SSH keys


2. clone the bitbucket repo :: 
    
    $ git clone git@bitbucket.org:py_gem/pyclusterevol.git

Now tell git who you are, and that you have an access to the repository. ::

    $ git init... I forgot! :)

Then to/from which branch you will push/pull. ::

    $ git --set-upstream remote_name branch_name


Revert recent modifications
---------------------------

Discard modification to a file/files ::

    $ git checkout -- file.name(s)

For example, if you made a change to files on a server and pushed the modification.
But forgotten doing so, you again start to modify a file on your laptop and realize that you are doing something redundant and you could simply pull from the repository. Then you can use "chekcout --" to revert the file so that you can pull from the repo.


Remove last commit ::

    $ git reset --hard HEAD~1
    
    or

    $ git reset --hard HEAD~10

to remove last 10 commits.

Overwrite local with remote
---------------------------
::

    $ git fetch --all
    $ git reset --hard remote/branch (i.e., origin/Hoseung)
    $ git pull 



Show modifications made to the current commit
---------------------------------------------
::

    $ git diff


Retrieve a file from previous commit
------------------------------------
::
	$ git checkout HEAD~n deleted-file.name

