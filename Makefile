NAME		= shallow_water
CC		= g++
SRC		= main.cpp Vector3.cpp
CFLAGS		= -Wall	-W -Werror

IFLAGS=	-I/usr/include -I/usr/include/GL -I/usr/X11R6/include -I/usr/X11R6/include/GL

LFLAGS=	-L/usr/lib -L/usr/X11R6/lib -lX11 -lGL -lGLU -lglut

OBJ   =	$(SRC:.cpp=.o)

all :		$(NAME)

$(NAME) :	$(OBJ)
		$(CC) $(OBJ) $(LFLAGS) -o $(NAME)

.cpp.o :
		$(CC) $(CFLAGS) $(IFLAGS) $< -c -o $@

clean :
		$(RM) $(OBJ)
		$(RM) $(NAME) 
