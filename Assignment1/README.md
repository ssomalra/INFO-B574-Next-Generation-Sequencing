# Assignment 1
This assignment familiarizes students with IU supercomputers and the use of the terminal.

## Assignment Submission
**_1. Log into Quartz using the terminal. Go into your slate directory. Echo "Hello World"_**
```
cd /N/slate/ssomalra/
echo "Hello World"
```
Output: ```Hello World```

**_2. Load python version 3.9.8 in Quartz and run python._**
```
module load python/3.9.8
python
```
Output:

<img width="600" alt="Screen Shot 2024-07-12 at 12 55 26 PM" src="https://github.com/user-attachments/assets/3ede1c61-3682-419c-a7c3-a7fc74f9364e">

<br> **_3. Generate a test file using the following command:_** ```for num in {a..z}; do printf $num"\tyourusername\n"; done | nl -n ln -w1 > test.txt```

Output:

<img width="300" alt="Screen Shot 2024-07-12 at 12 57 32 PM" src="https://github.com/user-attachments/assets/3620789d-5cca-4bfe-84ae-d0c119fe60d1">

<br> **_4. Add a fourth column concatenating the first two columns (The delimiter is "\t"). Write a script to generate the outcome using awk from the generated test.txt._**

```
awk -vOFS='\t' '{$4 = $1$2}1' test.txt > test2.txt
```

Output:

<img width="250" alt="Screen Shot 2024-07-12 at 1 01 20 PM" src="https://github.com/user-attachments/assets/9871d545-d4cc-448e-841a-932553f68b82">
