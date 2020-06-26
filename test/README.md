This folder contains subfolders with several types of tests:
+ *fast* contains verification tests that are used in Continuous Integration and should always pass on production branches
+ *bulk* contains example case that should be migrated towards another folder some day, an attempt is paid that the tests in bulk pass on production branches
+ *cases* contains work on real data (not versioned on the main git repository) the folder is as is and may not run with the current production version - check the git log of the files you are interested in
+ *bugs* should hopefuly be empty but is a transitory place to store scripts used to find a bug. It is a good practice to transform this scripts into generic test when the bug is found.
+ *linalg* contains tests that check the behaviour of LinearSolver features in real-test conditions and with stand-alone tests
