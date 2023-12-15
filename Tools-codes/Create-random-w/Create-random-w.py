def create_random(path,number,beginning,end):
    import random
    with open(path, "w") as file:
        for _ in range(number):
            random_value = random.uniform(beginning, end)
            file.write(str(random_value) + "\n")
    pass

def main():
    number_of_node=1000
    beginning_of_range=-0.5
    end_of_range=0.5
    address_name_for_save='./'+str(number_of_node)+','+str(beginning_of_range)+','+str(end_of_range)+'.txt'
    create_random(address_name_for_save,number_of_node,beginning_of_range,end_of_range)
    print("Done :)")
    pass


if __name__=="__main__":
    main()