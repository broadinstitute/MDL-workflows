version 1.0

# Temporary WDL 1.0 placeholder so Dockstore can parse and publish the entry.
# The real workflow (WDL 1.1) is restored in the next commit.
workflow sharding_for_CUDLL_scattered {
    input {
        String sample_name
    }

    call Placeholder { input: sample_name = sample_name }

    output {
        String message = Placeholder.message
    }
}

task Placeholder {
    input {
        String sample_name
    }

    command <<<
        echo "placeholder for ~{sample_name}"
    >>>

    output {
        String message = read_string(stdout())
    }

    runtime {
        docker: "ubuntu:22.04"
    }
}
