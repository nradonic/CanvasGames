<?php

if ($_SERVER['REQUEST_METHOD'] != "POST") {
    die("Improper request method");
}

if (!empty($_POST)) {
    foreach ($_POST as $key => $field) {
        if (strlen($field) === 0) {
            $missing_input[] = $key;
        }
    }
}   

if (!empty($missing_input)) {
    die("Error: " . implode(", ", $missing_input));
}




?>