# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    Makefile                                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: ashea <ashea@student.21-school.ru>         +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2021/07/04 11:30:46 by ashea             #+#    #+#              #
#    Updated: 2021/07/04 11:30:48 by ashea            ###   ########.fr        #
#                                                                              #
# **************************************************************************** #

CC = clang++
CFLAGS = -Werror -Wall -Wextra -I .
POLINOM = polinom
POLINOM_SRCS = polinom.o

%.o : %.c
	@$(CC) $(FLAGS) $< -c

all : $(POLINOM) 

$(POLINOM) : $(POLINOM_SRCS)
	@$(CC) $(POLINOM_SRCS) -o $(POLINOM)
	@printf "\e[32m$@ built\e[0m\n"

clean :
	@rm -f $(POLINOM_SRCS)
	@printf "\e[31mclean done\e[0m\n"

fclean : clean
	@rm -f $(POLINOM) 
	@printf "\e[31mfclean done\e[0m\n"

re : fclean all

.PHONY: all clean fclean re
