#include <cstdio>
#include <unistd.h>

int main(int argc, char* argv[]) {
	if(argc != 3) {
		printf("\x1b[31mUsage: ./run_interaction [grader] [user]\x1b[0m\n");
	}

	int fd_grader[2]; // grader -> user channel
	int fd_user[2]; // user -> grader channel

	if(pipe(fd_grader) != 0) printf("\x1b[31;1mError: Grader pipe failed\x1b[0m\n");
	if(pipe(fd_user) != 0) printf("\x1b[31;1mError: User pipe failed\x1b[0m\n");

	int childid;
	if((childid = fork()) == -1)
		printf("\x1b[31;1mError: Fork failed\x1b[0m\n");

	if(childid == 0) // grader
	{
		dup2(fd_grader[1], fileno(stdout));
		dup2(fd_user[0], fileno(stdin));
		close(fd_grader[1]);
		close(fd_user[0]);

		execl(argv[1], argv[1]);

		return 0;
	}

	// user
	dup2(fd_user[1], fileno(stdout));
	dup2(fd_grader[0], fileno(stdin));
	close(fd_user[1]);
	close(fd_grader[0]);

	execl(argv[2], argv[2]);
	return 0;
}
