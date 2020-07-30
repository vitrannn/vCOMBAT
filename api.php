<?php

define("EXECUTABLE_PATH", "bin/tuberculosis_simulation");
define("STD_OUT_FILE_PATH", "bin/stdout.out");
define("OUT_FILE_PATH", "bin/output.out");
@$useJSONP = $_GET['callback']; //checks if the call is being made with jsonp. if true, returns json output as a text within given callback function
$dictionary = false; //returns output as object dictionary. otherwise, double array without columns
//eg. [{'a':12,'at':123,'t':123,'l0':123,'l1':123},...]

function exception_handler($exception) {
    http_response_code(500);
    echo "Error: " . $exception->getMessage() . "\n";
}

function delete_out_files() {
    if (file_exists(STD_OUT_FILE_PATH)) {
        unlink(STD_OUT_FILE_PATH);
    }
    if (file_exists(OUT_FILE_PATH)) {
        unlink(OUT_FILE_PATH);
    }
}

set_exception_handler("exception_handler");

// 1. parse input parameters and construct exec command string:
$input = json_decode(file_get_contents('php://input'), true);

$options = array();
if ($input) {
    $allowedOptions = array( 'V', 'n', 'r', 'k', 'R', 'K', 'A', 'D', 'C', 't', 'd', 'p', 'S', );    // for more information on options refer to "tuberculosis_simulation" documentation
    foreach ($allowedOptions as $allowedOption) {
        if (array_key_exists($allowedOption, $input)) {
            $options[$allowedOption] = $input[$allowedOption];
        }
    }
}
$options['m'] = OUT_FILE_PATH;

$command = EXECUTABLE_PATH;
foreach ($options as $name => $value) {
    $command .= " -" . $name . " " . $value;
}
$command .= " -v > " . STD_OUT_FILE_PATH; // add verbose and redirect stdout to file

// 2. execute command:
//delete_out_files();    // delete old results, if any
exec($command);

// 3. parse results:
$responseData = array();

// parse verbose std output:
$stdOutContent = @file_get_contents(STD_OUT_FILE_PATH);
if ($stdOutContent === FALSE) {
    throw new Exception("Unable to open file: " . STD_OUT_FILE_PATH);
}
$responseData["verbose"] = $stdOutContent;

// parse "-m" file output:
$file = @fopen(OUT_FILE_PATH, "r");
if ($file === FALSE) {
    throw new Exception("Unable to open file: " . OUT_FILE_PATH);
}
$output = array();
if($dictionary){
    while(!feof($file)) {
        $str = fgets($file);
        $array = array();
        foreach(str_getcsv($str, " ") as $strValue) {
            if (!empty($strValue)) {
                array_push($array, floatval($strValue));
            }
        }
        if (!empty($array)) {
            array_push($output, $array);
        }
    }
}else{
    while(!feof($file)) {
        $str = fgets($file);
        $array = array();
        foreach(str_getcsv($str, " ") as $strValue) {
            if (!empty($strValue)) {
                array_push($array, floatval($strValue));
            }
        }
        if (!empty($array)) {
            array_push($output, $array);
        }
    }
}
$responseData["output"] = $output;
fclose($file);

//delete_out_files();

// 4. respond:
http_response_code(200);
if($useJSONP){ //if jsonp, return in given callback function
   echo "$useJSONP(".json_encode($responseData).')';
}else{ //return as regular json file
    header("Content-type:application/json; charset=utf-8");
    echo json_encode($responseData);
}

?>
