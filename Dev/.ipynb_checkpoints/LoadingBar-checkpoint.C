// C program to create loading bar
#include <stdio.h>
  
// Function to creating loading barvoid 
void loadingBar()
{

    char b = '#';
    cout << "\n\n\n\n";
    cout << "\n\n\n\n\t\t\t\t\t"
           << "Loading...\n\n";
    cout <<"\t\t\t\t\t";
  
    // Print initial loading bar
    for (int i = 0; i < 26; i++)
        cout << ("%c", b);
  
    // Set the cursor again starting
    // point of loading bar
    cout << "\r";
    cout << "\t\t\t\t\t";
    // Print loading bar progress
    for (int i = 0; i < 26; i++) {
        cout << ("%c", b);
  
        // Sleep for 1 second
    }
}
  
// Driver Code
int LoadingBar()
{
    // Function Call
    loadingBar();
    return 0;
}