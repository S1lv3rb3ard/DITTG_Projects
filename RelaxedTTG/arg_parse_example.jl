using ArgParse


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "theta12"
            help = "first twist angle theta12"
            required = true
        "theta23"
            help = "second twist angle theta23"
            required = true
end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    println(parsed_args["theta12"])
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
end

main()