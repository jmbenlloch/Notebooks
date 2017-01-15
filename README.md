# Notebooks
This is the invisible cities (IC) Notebook repository.

## Working with Notebooks. The IC recipe

Notebooks (Nbs)  are a powerfull tool for development. In practice, they can be used to draft and quickly develop the ideas that will eventually coalesce into reconstruction or analysis code. At the same time, Nbs can be seen as "live" documentation. Thus, the workings of an algorithm, reconstruction or analysis code can be very effectively described using the capabilities of Nbs to combine executable code (which in particular can produce images shown in the Nbs) and rich text.

To understand the workflow associated to Nbs, imagine that a developper needs to write a *complex algorithm* which requires several functions and eventually the manipulation of complex types. She creates the Nb complex_algo.ipynb to quickly draft those functions and classes, test the algorithm for particular cases, and try ideas that allow her to design some robust tests. At this stage the Nb does not need to look pretty, or be stable. Indeed, it is expected that the Nb changes continuously, as she tries different ideas and methods.

Therefore, complex_algorithm.ipynb is only useful for our heroine and consequenty lives only on her own fork of Notebooks. In fact, it is entirely up to her to write a single or many Nbs to study and solve her problem. She may choose, for example, to write a complex_algo.ipynb Nb and an a
complex_algo_tests.ipynb. There are no inherent restrictions, since the Nbs live on her fork, and in fact are modified only her.

It is important to understand the Nbs do not allow interactive collaboration. In practice, one cannot use git or any other similar tool to allow different users to work concurrently in the same Nb. This is due to the fact that the Nb will translate variables suh as the adress of a file or image, to the specific local data of each user, resulting in conflicts as soon as two users try to work in the same Nb. For this reason any user who is not the owner of the notebook must stash the changes before comitting.

For this very same reason, Nbs are not valid as representation of reconstruction or analysis code. Once the developer has understood the coding of *complex_algrithm* she writes a python (and or cython) script complex_algorithm.py *and* the corresponding test suite complex_algorithm_test.py and submits a pull request to nextic/IC. Unlike the Nbs, the python scripts are apt to collaborative work, and so a cycle of improvement may follow the initial commit, both at the level of the algorithm and of the associated tests and auxiliary functions.

Of course, she can feed back to her Nbs, the changes that occur as the pyhon (cython) scripts that code her algorithm are examined by other developers. Thus, the Nbs in her fork keep changing, although, eventually, they become more stable, as early draft code is removed, and code that lives in the python scripts is calledfrom the Nb rather than embeded in the Nb itself.

Eventually, the Nb may become an excellent description of how *difficult problem* was solved by
*complex algorithm* (ideas and notions can be coded in rich text), including live demonstration of performance. At this point, complex_algorithm.ipynb has became a *document* of general interest and she can decide to make a pull request to nextic/Notebooks. The central repository, then, is supposed to include only those clean, mature (concise) live documentation NBs.

## The structure of the repository

The Nb repository follows the same structure than the IC repository, with a manage.sh script that needs to be executed (to define the python version and the environment variable NBDIR which points to the top of the directory), then directories hanging from invisible_cities (cities, core, etc). There is one directory called files, where the user can place *demonstration* files (which are expected to be small, since they are eventually pushed to remote repositories). 

## Working with notebooks & git

1. Keep all your notebooks in your fork
2. Only consider uploading a notebook to the central repository when you have a mature version.
3. To upload to the central repository you must follow the usual git/IC convention summarize in: https://github.com/nextic/IC/blob/master/CONTRIBUTING.md
4. In short: create a topic branch, commit your candidate notebook to it, push it to your fork, make a pull request.
5. When/if the pull request is approved your flamant notebook is part of the central repository **and you are the owner of it**.
6. If you are the owner of a notebook you can modify it following the procedure we have just described.
7. If you are not the owner of a notebook be aware that:
7.1 Any attempt to modify the notebook will result in a conflict.
7.2 Just the fact of openning and running the notebook will modify it. Git will detect the modification and will ask if you want to stage it. **You can not**. At this stage the only option is to stash the notebook.
7.3 Alternatively (and probably a better strategy) make a copy of the notebook you just downloaded **without opening it**. To be specific, if you just got the notebook diomira.ipynb from the central repository and your user is braveheart then make a copy diomira_braveheart.ipynb. This new notebook can stay in your fork and you can do with it as you please. Importantly, since you have not modified the original diomira.ipynb you don't need to stash changes.

May the force be with you.