# Pre-processing module


def main():
    
    # Retrieve and merge junctions
    merged_junctions = process_and_merge_junctions(csv_path, design_parameters)

    # Stage reference genome if needed

    # For each junction, calculate the design region
    for j in merged_junctions:
        regions = calculate_primer_design_regions(j[0], j[1], design_parameters)



# Run main
if __name__=="__main__":
    main()