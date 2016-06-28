
Getting Started
============================================================================

Prerequisites
---------------------------------------------------------------------------
Some python packages are needed.


I recommend using virtualenv.::
   $ virtualenv venv
Where 'venv' is the name of a new directory inside which packages are installed.



As you can see in the above equation, with the *x* variable representing our genome list of integers, the *f(x)* shows our evaluation function, which is the sum of '0' values in the list. For example, if we have a list with 10 elements like this: ::
   
   x = [1, 2, 3, 8, 0, 2, 0, 4, 1, 0]


we will get the raw score [#rawscore]_ value of 3, or *f(x)* = 3. It's very simple to understand. Now, let's code this.

We will define our :term:`evaluation function` **"eval_func"** as: ::

   # This function is the evaluation function, we want
   # to give high score to more zero'ed chromosomes
   def eval_func(chromosome):
      score = 0.0

      # iterate over the chromosome elements (items)
      for value in chromosome:
         if value==0:
            score += 1.0
      
      return score

As you can see, this evaluation function tests each element in the list for equality with '0' and returns the proportional score value. The :class:`G1DList.G1DList` chromosome is not a python list by itself but it encapsulates one and exposes the methods for this list, like the iterator used in the above loop.
The next step is the creation of a :term:`sample genome` [#samplegenome]_ for the Genetic Algorithm. We can define our genome as this: ::

   # Genome instance
   genome = G1DList.G1DList(20)

   # The evaluator function (objective function)
