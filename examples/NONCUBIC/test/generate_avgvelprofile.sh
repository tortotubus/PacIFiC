#!/usr/bin/env bash

input="velocity.txt"
output="avgvelprofile.txt"

awk -v output="$output" '
BEGIN {
    line_count = 0
    total_time = 0.0
    header_line = ""
}

{
    # If this is the first line starting with "t" and header_line is not yet captured, store it
    if ($1 == "t" && header_line == "") {
        header_line = $0
    }

    # Determine time and velocity fields
    if ($1 == "t" && NF > 1) {
        # Line format: t <time> <vel1> <vel2> ...
        t = $2
        start_idx = 3
        vcount_current = NF - 2
    } else {
        # If not starting with "t", assume format: <time> <vel1> <vel2> ...
        if (NF > 1) {
            t = $1
            start_idx = 2
            vcount_current = NF - 1
        } else {
            # Not a valid data line
            next
        }
    }

    # Process only lines where time > 50
    if (t > 50) {
        if (line_count == 0) {
            # First valid line sets the number of velocities and initializes sums
            vel_count = vcount_current
            for (i = 1; i <= vel_count; i++) {
                sum[i] = 0.0
            }
        } else {
            # Ensure consistent number of velocities in subsequent lines
            if (vcount_current != vel_count) {
                print "Inconsistent number of velocities detected." > "/dev/stderr"
                exit 4
            }
        }

        total_time += t

        # Accumulate velocity sums
        for (i = 0; i < vel_count; i++) {
            sum[i+1] += $(start_idx + i)
        }

        line_count++
    }
}

END {
    if (line_count == 0) {
        print "No lines found with t > 50." > "/dev/stderr"
        exit 2
    }

    # Compute averages
    avg_time = total_time / line_count

    # Write the first line that started with "t" to the output (as a header)
    print header_line > output

    # Append the averaged values to the output
    printf("%.6f", avg_time) >> output
    for (i = 1; i <= vel_count; i++) {
        printf(" %.6f", sum[i]/line_count) >> output
    }
    printf("\n") >> output
}' "$input"
