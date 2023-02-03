module NARDS

export nards, plot_results

using Statistics, ProgressMeter, LaTeXStrings, FASTX, Random, ThreadsX


"""
    const DEFAULT_CATEGORIES = Dict{String, Array{Char}}(
        "mu" => ['S' 'T' 'N' 'Q' 'C' 'H'], 
        "h" => ['I' 'L' 'M' 'V'], 
        "+" => ['R' 'K'], 
        "-" => ['E' 'D'],
        "pi" => ['F' 'W' 'Y'], 
        "A" => ['A'], 
        "P" => ['P'], 
        "G" => ['G'], 
    )

The default way to categorise the amino acids for the analysis of the clusters. 
"""
const DEFAULT_CATEGORIES = Dict{String, Array{Char}}(
    "mu" => ['S' 'T' 'N' 'Q' 'C' 'H'], 
    "h" => ['I' 'L' 'M' 'V'], 
    "+" => ['R' 'K'], 
    "-" => ['E' 'D'],
    "pi" => ['F' 'W' 'Y'], 
    "A" => ['A'], 
    "P" => ['P'], 
    "G" => ['G'], 
)

"""
    const DEFAULT_CATEGORY_ORDER = Array{String}(["mu", "h", "+", "-", "pi", "A", "P", "G"])

The default order of the categories in the resulting heat map. 
"""
const DEFAULT_CATEGORY_ORDER = Array{String}(["mu", "h", "+", "-", "pi", "A", "P", "G"])

"""
    const DEFAULT_CATEGORY_LABELS = Array{LaTeXString}([L"\\mu", "h", "+", "-", L"\\pi", "A", "P", "G"])

The default labels for the categories in the heatmap output. 
"""
const DEFAULT_CATEGORY_LABELS = Array{LaTeXString}([L"\mu", "h", "+", "-", L"\pi", "A", "P", "G"])

"""
    nards(
        seq::Array{Char}; 
        num_null_models::Int=100_000, 
        g::Int=5, 
        categories::Dict{String, Array{Char}}=DEFAULT_CATEGORIES, 
        category_order::Array{String}=DEFAULT_CATEGORY_ORDER, 
        show_progress_bar::Bool=false,
    )

Performs an analysis on `seq` using the NARDINI algorithm. 

There are several key differences 

# Arguments:
- `seq`: The sequence of amino acids that will be analyzed. 
- `num_null_models`: The number of null models to use for computing the z-scores.
- `g`: The size of the sliding window in the asymmetry calculation. 
- `categories`: A `Dict` that contains the categories we will group the amino acids into.
- `category_order`: The order the categories will be displayed in the heat map. (Must contain all keys in `categories`.)
- `show_progress_bar`: Whether or not to display a progress bar on the command line. 
"""
function nards(
        seq::Array{Char}; 
        num_null_models::Int=100_000, 
        g::Int=5, 
        categories::Dict{String, Array{Char}}=DEFAULT_CATEGORIES, 
        category_order::Array{String}=DEFAULT_CATEGORY_ORDER, 
        show_progress_bar::Bool=false
    )

    random_shuffles = fill('-', (num_null_models, length(seq)));
    Threads.@threads for i in 1:num_null_models
        random_shuffles[i, :] = Random.shuffle(seq); 
    end

    kappas = zeros(8,8)

    # Iterate over the category combinations. 
    total_categories = 0
    for (i, l1) in enumerate(category_order)
        for (j, l2) in collect(enumerate(category_order))[i:length(categories)]
            total_categories += 1
        end
    end
    if show_progress_bar
        prog = Progress(total_categories)
    end

    for (i, l1) in enumerate(category_order)
        for (j, l2) in collect(enumerate(category_order))[i:length(categories)]
            # Compute the categories. 
            if i == j   # Omega type 
                cat1 = categories[l1]
                cat2 = Array{Char}(undef, 0)
                for l in keys(categories) 
                    if l != l1 
                        append!(cat2, categories[l]) 
                    end
                end
            else        # Delta type
                cat1 = categories[l1]
                cat2 = categories[l2]
            end

            # Skip this one (set to zero) if there is less than 10% composition. 
            if ThreadsX.sum(ThreadsX.map((x) -> x in cat1, seq)) / length(seq) < 0.1 || ThreadsX.sum(ThreadsX.map((x) -> x in cat2, seq)) / length(seq) < 0.1
                kappas[i,j] = 0.0
                if show_progress_bar
                    next!(prog) # Update the progress bar. 
                end
                continue
            end
            
            # The asymmetry function. 
            function delta(seq)
                # Compute the asymmetry. 
                function asymmetry(seq)
                    f1 = sum(map((x) -> x in cat1, seq)) / length(seq)
                    f2 = sum(map((x) -> x in cat2, seq)) / length(seq)
                    # println(map((x) -> x in cat1, seq))
                    if f1 + f2 == 0 
                        return 0.0
                    else 
                        return (f1 - f2)^2 / (f1 + f2)
                    end
                end
                global_asymmetry = asymmetry(seq)
                local_asymmetries = ThreadsX.map(asymmetry, ((seq[i:i+g-1]) for i in 1:length(seq)-g+1))
                ThreadsX.map!((x) -> (x - global_asymmetry)^2, local_asymmetries, local_asymmetries)
                return mean(local_asymmetries)
            end

            kappa_values = ThreadsX.map(delta, eachrow(random_shuffles))
            # println(l1, '\t', l2,'\t', deltas[1:10]) # Print the first ten deltas for debugging. 

            function f(x) 
                if x in cat1
                    return -1
                elseif x in cat2
                    return 1
                else 
                    return 0
                end
            end
            max_delta = delta(sort(seq, by=f))

            ThreadsX.map!((x) -> x / max_delta, kappa_values, kappa_values)

            m = mean(kappa_values)
            s = std(kappa_values)
            q = delta(seq)/max_delta

            kappas[i,j] = (q-m)/s # Save the z-score. 

            if show_progress_bar
                next!(prog) # Update the progress bar. 
            end
        end
    end

    replace!(kappas, NaN=>0.0) # Replace NaNs from things that have zero options. 

    return kappas
end


using CairoMakie
function plot_results(z_scores::Array{Float64, 2}, path::String; category_labels::Array{LaTeXString}=DEFAULT_CATEGORY_LABELS)
    fig = Figure(resolution=(800, 600))
    ax = Axis(fig[1, 1],
        xticks = (1:length(category_labels), category_labels), 
        yticks = (1:length(category_labels), category_labels),
        yreversed = true
    )
    hm = heatmap!(fig[1,1], transpose(z_scores), colormap=:bwr, colorrange=(-3, 3))
    cb = Colorbar(fig[1, 2], hm)

    # Sluging for safety. 
    path = replace(path, '/' => '_'); 
    path = replace(path, '\\' => '_'); 

    save(path, fig)
end


#= 
    Below is the CLI/IO component of the application. 
=#
using ArgParse, DataFrames, CSV

function analyze_seq(seq::AbstractString; g::Int64=5, output_path::AbstractString=".", num_null_models::Int=100_000, title::AbstractString="", show_progress_bar::Bool=false)
    if show_progress_bar
        println("Analyzing sequence: ", seq)
        println("\tWindow width:", g)
        println("\tSaving output to:", output_path)
    end

    # Generate the z-score matrix. 
    z_scores = nards(collect(seq), num_null_models=num_null_models, g=g, show_progress_bar=show_progress_bar)

    # Generate the z-score heatmap. 
    safe_title = replace(title, '/' => '_')
    if title == "" 
        plot_results(z_scores, string(output_path, "/heatmap.png")) # Default file name. 
    else 
        plot_results(z_scores, string(output_path, '/', safe_title, ".png")) # Chosen file name. 
    end

    # Generate the z-score CSV. 
    # TODO: Generate the CSV. How to do this???? Probably with a DF. 
    df = DataFrame(z_scores, :auto)

    rename!(df,:x1 => :mu)  # Hard code this for now. 
    rename!(df,:x2 => :h)   # TODO: make this not hard coded lol. 
    rename!(df,:x3 => :+)
    rename!(df,:x4 => :-)
    rename!(df,:x5 => :pi)
    rename!(df,:x6 => :A)
    rename!(df,:x7 => :P)
    rename!(df,:x8 => :G)
    
    # display(df)
    if safe_title == "" 
        p = string(output_path, "/z_scores.csv") # Default file name. 
    else 
        p = string(output_path, '/', safe_title,  "_z_scores.csv") # Chosen file name. 
    end
    CSV.write(p, df)

    return z_scores
end

function parse_cli()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "input"
            help = "The sequence (or fasta) to be analyzed"
            required = true
        "--window_size", "-g"
            help = "The size of the sliding window to use"
            arg_type = Int
            default = 5
        "--outout_path", "-o"
            help = "The path to the file where the output will be stored"
            arg_type = String
            default = "."
        "--num_null_models", "-n"
            help = "The number of sequence scrambles to use in the null model"
            arg_type = Int
            default = 100_000 
        "--title" 
            help = "The title of the output plot"
            arg_type = String
            default = "" 
        "--seed", "-s"
            help = "The random seed that will be used to generate the null models for reproducability. If set to zero then one will be selected randomly."
            arg_type = Int
            default = 0
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_cli()

    # Check to see if we have set a random seed. 
    if parsed_args["seed"] != 0 
        Random.seed!(parsed_args["seed"])
    end
    
    if isfile(parsed_args["input"]) # We are processing a fasta. 
        input = parsed_args["input"]
        g = parsed_args["window_size"]
        output_path = parsed_args["outout_path"]
        num_null_models = parsed_args["num_null_models"]
        title = parsed_args["title"]

        # Reader for each of the sequences.
        f = open(input) 
        reader = FASTAReader(f)
        seqs = collect(reader)
        close(reader)

        # Record all of the z-scores. 
        z_scores = zeros(length(DEFAULT_CATEGORY_LABELS), length(DEFAULT_CATEGORY_LABELS), length(seqs))

        # Process each of the sequences in the fasta. 
        prog = Progress(length(seqs), desc="Processing FASTA input")
        Threads.@threads for (i, entry) in collect(enumerate(seqs))
            s = sequence(entry)
            n = description(entry)

            zs = analyze_seq(s, 
                    g=g, 
                    output_path=output_path, 
                    num_null_models=num_null_models, 
                    title=n, 
                    show_progress_bar=false
                )

            z_scores[:, :, i] = zs
            next!(prog) # Update the progress bar. 
        end

        # Generate the matrix and labels for the heat map. 
        total_categories = 0
        labels = LaTeXString[]
        category_order = DEFAULT_CATEGORY_ORDER # Hard-coded for now. 
        categories = DEFAULT_CATEGORIES
        for (i, l1) in enumerate(category_order)
            for (j, l2) in collect(enumerate(category_order))[i:length(categories)]
                total_categories += 1
                l = DEFAULT_CATEGORY_LABELS[i] * L"-" * DEFAULT_CATEGORY_LABELS[j]
                push!(labels, l)
            end
        end 
        collective_matrix = zeros(length(seqs), total_categories)
        for k in range(1, length(seqs))
            cur_idx = 1
            for (i, l1) in enumerate(category_order)
                for (j, l2) in collect(enumerate(category_order))[i:length(categories)]
                    collective_matrix[k, cur_idx] = z_scores[i,j,k]
                    cur_idx += 1
                end
            end
        end 

        # Create and save the heatmap. 
        fig = Figure(resolution=(800, 1200))
        x_names = String[]
        for s in seqs
            push!(x_names, description(s))
        end
        ax = Axis(fig[1, 1],
            xticks = (1:length(x_names), x_names), 
            yticks = (1:length(labels), labels),
            yreversed = true, 
            xticklabelrotation=45.0
        )
        hm = heatmap!(fig[1,1], collective_matrix, colormap=:bwr, colorrange=(-3, 3))
        cb = Colorbar(fig[1, 2], hm)

        save(string(output_path, "/meta_heatmap.png"), fig)

    else # We are processing a sequence. 
        input = parsed_args["input"]
        g = parsed_args["window_size"]
        output_path = parsed_args["outout_path"]
        num_null_models = parsed_args["num_null_models"]
        title = parsed_args["title"]

        analyze_seq(input, 
                g=g, 
                output_path=output_path, 
                num_null_models=num_null_models, 
                title=title
            )
    end
end

# If we are running the script from the command line. 
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module NARDS
