API interface endpoint to "tuberculosis_simulation"
========================================


Requirements:
========================================
* PHP 5
* "tuberculosis_simulation" executable
* "inputsample.txt" file with data (one float value per line)


Setup:
========================================
* Organize the files as follows:
./
    bin/
        tuberculosis_simulation
    api.php
    inputsample.txt

* For information on compiling and linking of "tuberculosis_simulation" source code refer to "tuberculosis_simulation" documentation. Ideally the steps are:
$ cd project_root
$ cmake .
$ make
# then copy produced executable file


Input:
========================================
POST JSON request. Optional parameters (case sensitive):
{
    'V': intracellularVolume,
    'n': targetMoleculeCount,
    'r': replicationThreshold,
    'k': killingThreshold,
    'R': baselineReplicationRate,
    'K': maximumKillingRate,
    'A': targetAssociationRate,
    'D': targetDissociationRate,
    'C': carryingCapacity,
    't': time,
    'd': startingAntibiotic,
    'M': AntibioticMolWeight,
    'p': startingPopulation,
    'S': steppingFunction,
}

For more information refer to "tuberculosis_simulation" documentation.


Output:
========================================
On Success:
    status: 200
    body:
        JSON:
            {
                "verbose": string, the std console output "tuberculosis_simulation" produces; for more information refer to "tuberculosis_simulation" documentation.
                "output": 2-dimensional array of floats, the output file "tuberculosis_simulation" produces when called with "-m" key; for more information refer to "tuberculosis_simulation" documentation.
            }
On Failure:
    status: 500
    body:
        error message string


Run:
========================================
php -S localhost:8888 api.php


Test:
========================================
curl -X POST -H "Content-Type: application/json" -d '{ "d": 10, "p": 1000 }' "http://localhost:8888/"
