import numpy as np
import numpy.linalg as la

class Joint:
    def __init__(self, joint_id, loc_x, loc_y):
        self.joint_id = joint_id
        self.loc_x = loc_x
        self.loc_y = loc_y

    def __str__(self):
        return str(self.joint_id) + " (" + str(self.loc_x) + ", " + str(self.loc_y) + ")"

    def init_displacements(self,disp_x, disp_y):
        if disp_x:
            self.disp_x = 0.0
        if disp_y:
            self.disp_y = 0.0

    def init_reactions(self,react_x, react_y):
        if react_x:
            self.react_x = 0.0
        if react_y:
            self.react_y = 0.0

    def init_force(self, force_x, force_y):
        self.force_x = force_x
        self.force_y = force_y

class Reaction:
    def __init__(self, joint, react_x, react_y):
        self.joint = joint
        if react_x:
            self.react_x = 0.0
        if react_y:
            self.react_y = 0.0

    def __str__(self):
        return self.joint.__str__() + " x fixed: " \
        + str(self.react_x) + " y fixed: " + str(self.react_y)

class Force:
    def __init__(self, joint, force_x, force_y):
        self.joint = joint
        self.force_x = force_x
        self.force_y = force_y

    def __str__(self):
        return self.joint.__str__() + " x force: " \
        + str(self.force_x) + " y force: " + str(self.force_y)

class Disp:
    def __init__(self,joint, fixed_x, fixed_y):
        self.joint = joint
        if not fixed_x:
            self.disp_x = 0.0
        if not fixed_y:
            self.disp_y = 0.0



class Member:
    def __init__(self, member_id, joint_i, joint_j, area, E):
        self.member_id = member_id
        self.joint_i = joint_i
        self.joint_j = joint_j
        self.area = area
        self.E = E
        self.length = ((self.joint_j.loc_x - self.joint_i.loc_x) ** 2 + (self.joint_j.loc_y - self.joint_i.loc_y) ** 2) ** .5
        self.alpha = (self.joint_j.loc_x - self.joint_i.loc_x) / self.length
        self.beta = (self.joint_j.loc_y - self.joint_i.loc_y) / self.length
        self.K = np.array((self.alpha**2,self.alpha*self.beta,self.alpha*self.beta,self.beta**2)).reshape(2,2) \
         * self.area * self.E / self.length

    def __str__(self):
        return str(self.member_id) + " " + str(self.joint_i.joint_id) + " " +  \
               str(self.joint_j.joint_id) + " " +  str(self.length) + " " +  str(self.alpha) + " " +  str(self.beta)

    def calc_member_forces(self):
        pass


class Global_Stiffness_Matrix:
    def __init__(self,joints):
        self.gsm = np.zeros((len(joints) * 2,len(joints) * 2))

    def __add__(self,member):
        ix = member.joint_i.joint_id * 2
        iy = member.joint_i.joint_id * 2 + 1
        jx = member.joint_j.joint_id * 2
        jy = member.joint_j.joint_id * 2 + 1
        ind = [ix,iy,jx,jy]
        member_k = np.hstack((member.K, - member.K))
        member_k = np.vstack((member_k, np.hstack((-member.K, member.K))))
        for i in range(len(ind)):
            for j in range(len(ind)):
                #print member_k[i][j]
                #print ind[i], ind[j]
                self.gsm[ind[i],ind[j]] += member_k[i][j]
                #print member.member_id
                #print self.gsm
        #print member_k
        #print ind
    def __str__(self):
        return str(self.gsm)

def init_joints():
    print "Joints**********"
    jointsA = np.zeros(())
    jointsA =  np.vstack((jointsA, np.array([0.0, 0.0])))
    jointsA =  np.vstack((jointsA, np.array([7.5,0.0])))
    jointsA =  np.vstack((jointsA, np.array([4.8,3.6])))
    print jointsA[1:,]
    joints = []
    joints.append(Joint(0,0.0,0.0))
    joints.append(Joint(1,7.5,0.0))
    joints.append(Joint(2,4.8,3.6))

    for i in joints:
        print i
    return joints

def init_forces():
    forces = []
    forces.append(Force(2,120,80))
    print "Forcess**********"
    for i in forces:
        print i
    return forces


def init_reactions():
    reactions = []
    reactions.append(Reaction(joints[0],True,True))
    reactions.append(Reaction(joints[1],True,True))
    print "Reactions**********"
    for i in reactions:
        print i
    return reactions

def init_members():
    members = []
    members.append(Member(0,joints[0], joints[2], 900, 200))
    members.append(Member(1,joints[2], joints[1], 900, 200))
    print "Members**********"
    for i in members:
        print i
        print i.K
    return members

def init_Global_Stiffness_Matrix(members):
    print "Global_Stiffness_Matrix **********"
    GK = Global_Stiffness_Matrix(joints)
    GK + members[0]
    GK + members[1]
    print GK
    return GK

if __name__ == '__main__':
    
    N = 3
    M = 2
    joints = np.zeros((N, 2))
    members = np.zeros((M,1))
    K = np.zeros((N*2, N*2))
    forces = np.zeros((N,2))
    joints = init_joints()
    forces = init_forces()
    reactions = init_reactions()
    members = init_members()
    GK = init_Global_Stiffness_Matrix(members)
    k = GK.gsm[4:6,4:6]
    f = np.array((forces[0].force_x, forces[0].force_y)).reshape(2,1)
    print f
    print k
    u = la.solve(k,f)
    #print la.inv(k)
    print u
    k=GK.gsm[0:4,4:6]
    print k
    print np.dot(k,u)




