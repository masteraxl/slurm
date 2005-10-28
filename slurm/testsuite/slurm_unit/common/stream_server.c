/* Slurm stream server test, read stdin until "quit" is entered, anything else is sent 
 * out to the specified SLURM port, anything received on that port is printed */
#include <src/common/log.h>
#include <netinet/in.h>
#include <src/common/slurm_protocol_api.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

extern int errno;
void* read2stdout_thread( void* arg )
{
	unsigned int buffer_len = 1024*1024 ;
	int bytes_read ;
	char buffer[ buffer_len ] ;
	slurm_fd *fd = (slurm_fd *) arg ;

	while(1)
	{
		if ( ( bytes_read = slurm_read_stream( *fd, buffer, buffer_len ) ) < 0 )
			break;
		buffer[ bytes_read ] = '\0' ;
		printf( "%s", buffer );
	}
	return NULL;
}

void stdout2socket_loop( slurm_fd fd )
{
	unsigned int buffer_len = 1024*1024 ;
	int bytes_read ;
	char buffer[ buffer_len ] ;
	char* current = buffer;
	int curr_size = 0 ;
	while (1)
	{
		if( (bytes_read = read(STDIN_FILENO, current, buffer_len-curr_size ) ) < 0 )
		{
			printf("read error\n");
			break;
		}
		
		curr_size += bytes_read ;
		if ( '\n' == buffer[ curr_size - 1 ] )
		{
			if ( strncmp( "quit", buffer, 4 ) == 0 )
				break;
			buffer[ curr_size ] = '\0' ;
			if ( slurm_write_stream( fd, buffer, curr_size ) < 0 )  /* fails */
			{
				printf("Could not send\n");
				break;
			}
			current = buffer ;
			curr_size = 0 ;
		}

	}
}

int main ( int argc , char * argv[] )
{
	log_options_t log_opts = { 1, LOG_LEVEL_DEBUG3,  LOG_LEVEL_INFO, LOG_LEVEL_QUIET } ;
	slurm_fd listen_socket ;
	slurm_fd worker_socket ;

	/* declare address structures */
	slurm_addr listen_address ;
	slurm_addr worker_address ;

	/* read thread */
	pthread_t read_pth ;
	int16_t port = 0;

	if (argc > 1)
		port = atoi( argv[1] ) ;

	if ((argc < 2) || (port < 1)) {
		printf("Usage: %s <port_number>\n", argv[0] );
		exit( 1 );
	}

	
	/* init address sturctures */
	log_init(argv[0], log_opts, SYSLOG_FACILITY_DAEMON, NULL);
	slurm_set_addr_uint ( & listen_address , port, SLURM_INADDR_ANY ) ;

	/* open and listen on socket */
	listen_socket = slurm_listen_stream ( & listen_address ) ;

	/* accept socket */
	worker_socket = slurm_accept_stream ( listen_socket , & worker_address ) ;

	if ( pthread_create( &read_pth, NULL, &read2stdout_thread, (void*) &worker_socket ) ) 	
	{
		printf("Could not create read_thread: error=%d\n", errno );
		exit( errno );
	}

	stdout2socket_loop( worker_socket );
	
	slurm_close_stream ( worker_socket ) ;
	slurm_close_stream ( listen_socket ) ;

	return 0 ;
}
