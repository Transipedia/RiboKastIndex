directory=$1
inputFile=$2
outputFile=$3

file_list=$(ls -p "$directory" | grep -v /)

header=("tag")

# Loop through each file name and add the modified filename (without extension) to the header array
for file_name in $file_list; do
    file_name=$(basename "$file_name")  # Extract only the filename
    file_name="${file_name%.*}"  # Remove the extension

    # Check if the filename is already in the header array
    if [[ ! " ${header[@]} " =~ " $file_name " ]]; then
        header+=("$file_name")
    fi
done

header_line=$(printf '%s\t' "${header[@]}")
header_line=${header_line%$'\t'}  # Remove the trailing tab

echo -e "$header_line" > header.tsv
cat header.tsv "$inputFile" > "$outputFile"
rm header.tsv
