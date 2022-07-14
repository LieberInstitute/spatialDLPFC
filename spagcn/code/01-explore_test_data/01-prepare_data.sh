#  Document how to retrieve SpaGCN's test data, and how to organize it inside
#  this repo.

#  Get absolute path to 'spython' repo
base_dir=$(git rev-parse --show-toplevel)

#  Manually download test data (done here on 9/22/21) from:
#  https://drive.google.com/drive/folders/1zten54vkjorp26T4iD0ApQGa9ut5eY42?usp=sharing
#  Download both folders, each of which are downloaded as zip files.

#  Upload each zip file to 'spagcn/raw-data/01-explore_test_data/'

#  Ignore the test data, which is too large to commit and push
echo "spagcn/raw-data/01-explore_test_data" >> $base_dir/.gitignore

#  Unzip the data
cd $base_dir/spagcn/raw-data/01-explore_test_data
unzip *.zip
rm *.zip
